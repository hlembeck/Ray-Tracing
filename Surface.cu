#include "Surface.cuh"

Face::Face(Point p1, Point p2, Point p3) : _p1(p1), _p2(p2), _p3(p3), _spectrum() {
	_orientation = Vector(p2 - p1).cross(Vector(p3 - p2));
}

bool Face::onPlane(Point p) {
	Vector v1(_p2 - _p1);
	Vector v2(_p3 - _p1);
	Vector v(p - _p1);
	return (v.dot(v1.cross(v2)) ? false : true);
}

bool Face::inTriangle(Point p) {
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
		return 0.0;
	t = (Vector(_p1 - p).dot(_orientation)) / t;
	return (inTriangle(p.translate(v, t)) ? t : 0.0);
}

__host__ Face* getStandardPrism() {
	Face* faces;
	cudaMallocManaged((void**)&faces, sizeof(Face) * 8);
	Point v1(-1.0, 0.0, 1.0);
	Point v2(-1.0, 1.0, 2.0);
	Point v3(-1.0, 0.0, 3.0);
	Point v4(1.0, 0.0, 1.0);
	Point v5(1.0, 1.0, 2.0);
	Point v6(1.0, 0.0, 3.0);
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

__device__ double Mesh::getIntersection(Ray r) {
	double intersection;
	for (unsigned int i = 0; i < _numFaces; i++) {
		if ((intersection = _faces[i].getIntersection(r.direction,r.origin)) > 0)
			return intersection;
	}
	return -1.0;
}