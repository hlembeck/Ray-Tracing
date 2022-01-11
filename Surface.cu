#include "Surface.cuh"

Face::Face(Point p1, Point p2, Point p3) : _p1(p1), _p2(p2), _p3(p3), _spectrum(PiecewiseSpectrum()) {
	_orientation = Vector(p2 - p1).cross(Vector(p3 - p2));
	_orientation.normalize();
}

Face::Face(Point p1, Point p2, Point p3, PiecewiseSpectrum s) : _p1(p1), _p2(p2), _p3(p3), _spectrum(s) {
	_orientation = Vector(p2 - p1).cross(Vector(p3 - p2));
	_orientation.normalize();
}

__device__ bool Face::onPlane(Point p) {
	Vector v1(_p2 - _p1);
	Vector v2(_p3 - _p1);
	Vector v(p - _p1);
	return (v.dot(v1.cross(v2)) ? false : true);
}

__device__ bool Face::inTriangle(Point p) {
	if (!onPlane(p))
		return false;
	Vector v1(_p1 - p);
	Vector v2(_p2 - p);
	Vector v3(_p3 - p);
	Vector t1 = v1.cross(v2);
	Vector t2 = v2.cross(v3);
	Vector t3 = v3.cross(v1);
	if (t1.dot(t2) < 0.0 || t2.dot(t3) < 0.0)
		return false;
	return true;
}

__device__ double Face::getIntersection(Vector v, Point p) {
	double t = v.dot(_orientation);
	if (t == 0.0)
		return -1.0;
	t = (Vector(_p1 - p).dot(_orientation)) / t;
	return (inTriangle(p.translate(v, t)) ? t : -1.0);
}

__device__ void Face::transform(double* m) {
	printf("in Face::transform \n");

	_p1.x = m[0] * _p1.x + m[1] * _p1.y + m[2] * _p1.z + m[3];
	_p1.y = m[4] * _p1.x + m[5] * _p1.y + m[6] * _p1.z + m[7];
	_p1.z = m[8] * _p1.x + m[9] * _p1.y + m[10] * _p1.z + m[11];

	_p2.x = m[0] * _p2.x + m[1] * _p2.y + m[2] * _p2.z + m[3];
	_p2.y = m[4] * _p2.x + m[5] * _p2.y + m[6] * _p2.z + m[7];
	_p2.z = m[8] * _p2.x + m[9] * _p2.y + m[10] * _p2.z + m[11];

	_p3.x = m[0] * _p3.x + m[1] * _p3.y + m[2] * _p3.z + m[3];
	_p3.y = m[4] * _p3.x + m[5] * _p3.y + m[6] * _p3.z + m[7];
	_p3.z = m[8] * _p3.x + m[9] * _p3.y + m[10] * _p3.z + m[11];
}

Face* getStandardPrism() {
	Face* faces;
	cudaMallocManaged((void**)&faces, sizeof(Face) * 8);
	Point v1(-1.0, 0.0, -1.0);
	Point v2(-1.0, 1.0, 0.0);
	Point v3(-1.0, 0.0, 1.0);
	Point v4(1.0, 0.0, -1.0);
	Point v5(1.0, 1.0, 0.0);
	Point v6(1.0, 0.0, 1.0);
	faces[0] = Face(v1, v2, v3);
	faces[1] = Face(v6, v5, v4);
	faces[2] = Face(v4, v5, v2);
	faces[3] = Face(v4, v2, v1);
	faces[4] = Face(v6, v4, v1);
	faces[5] = Face(v6, v1, v3);
	faces[6] = Face(v3, v2, v5);
	faces[7] = Face(v3, v5, v6);
	return faces;
}

Face* getStandardPlate(PiecewiseSpectrum& const s) {
	Face* faces;
	cudaMallocManaged((void**)&faces, sizeof(Face) * 2);
	Point v1(1.0,-1.0,0.0);
	Point v2(-1.0, 1.0, 0.0);
	Point v3(1.0, 1.0, 0.0);
	Point v4(-1.0, -1.0, 0.0);
	faces[0] = Face(v1,v2,v3,s);
	faces[1] = Face(v1,v4,v2,s);
	return faces;
}

void Face::print() {
	_p1.print();
	printf(" -> ");
	_p2.print();
	printf(" -> ");
	_p3.print();
	printf("\n");
	printf("Orientation: ");
	_orientation.print();
	printf("\n");
}

void Mesh::print() {
	for (unsigned int i = 0; i < _numFaces; i++) {
		_faces[i].print();
	}
}

__global__ void meshTransform(Face* faces, double* mat, unsigned int const len) {
	unsigned int x = threadIdx.x;
	unsigned int y = blockIdx.x;
	unsigned int i = y * 256 + x;
	if (i < len) {
		faces[i].transform(mat);
	}
}

void Mesh::transform(Transformation& const m) {
	unsigned int n = (_numFaces & 255 ? (_numFaces >> 8) + 1 : _numFaces >> 8);
	meshTransform<<<n, 256>>>(_faces, m.getTransformation(), _numFaces);
	cudaDeviceSynchronize();
}

__device__ void printArrayDEVICE(double* const arr, unsigned int const len) {
	for (unsigned int i = 0; i < len - 1; i++) {
		printf("%f, ", arr[i]);
	}
	printf("%f\n", arr[len - 1]);
};

__device__ bool Mesh::handleIntersection(Ray& r, double* s, double& angle) {
	//printf("in handleIntersection\n");
	double intersection, temp = DBL_MAX;
	for (unsigned int i = 0; i < _numFaces; i++) {
		intersection = _faces[i].getIntersection(r.direction, r.origin);
		if (intersection > 1.0 && intersection < temp) {
			//printf("%f\n", intersection);
			s = _faces[i]._spectrum.getPower();
			//printArrayDEVICE(s, STEPS);
			angle = r.direction.dot(_faces[i]._orientation);
			r.origin = r.origin.translate(r.direction,intersection);
			r.direction = r.direction - _faces[i]._orientation.scale(2 * (r.direction.dot(_faces[i]._orientation)));
			temp = intersection;
			//r.origin.print();
		}
	}
	return false;
}

__global__ void traceRay(Ray* rays, Mesh* meshes, unsigned int n, double** spectrums, double* angles) {
	int x = blockIdx.x;
	int y = blockIdx.y;
	for (unsigned int i = 0; i < n; i++) {
		if (meshes[i].handleIntersection(rays[x * raysPerPixel + y], spectrums[x * raysPerPixel + y], angles[x * raysPerPixel + y]))
			return;
	}
	angles[x * raysPerPixel + y] = 1.0;
}

void Scene::transform(Transformation& const t) {
	for (unsigned int i = 0; i < _numMeshes; i++) {
		_meshes[i].transform(t);
	}
}