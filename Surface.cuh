#include "LightColor.cuh"

const RefractiveIndex ambient{1,0};

class SurfaceProperties {
public:
	SurfaceProperties(double r, double r1, double r2) : _reflectance(r), _in(RefractiveIndex{ r1, r2}), _out(ambient) {};
	SurfaceProperties() : _reflectance(1.0), _in(RefractiveIndex{ 1,0 }), _out(ambient) {};
	~SurfaceProperties() {};
	double getInnerIndex(double l) { return _in.r1 + _in.r2 / (l * l); };
	double getOuterIndex(double l) { return _out.r1 + _out.r2 / (l * l); };
	double getReflectance() { return _reflectance; };

private:
	double _reflectance;
	RefractiveIndex _in;
	RefractiveIndex _out;
};

//Supports only a constant spectrum across the sphere.
class Sphere {
public:
	Sphere(Point o, double r, PiecewiseSpectrum s) : _origin(o), _radius(r), _spectrum(s) {};
	~Sphere() {};
private:
	Point _origin;
	double _radius;
	PiecewiseSpectrum _spectrum;
};

class Face {
	friend class Mesh;
public:
	Face(Point p1, Point p2, Point p3);
	Face(Point p1, Point p2, Point p3, PiecewiseSpectrum s);
	~Face() {};

	__device__ void transform(double* m);

	__device__ bool onPlane(Point p);
	__device__ bool inTriangle(Point p);
	//Returns t s.t p+vt lies on face, otherwise returns -1.0
	__device__ double getIntersection(Vector v, Point p);
	void print();
private:
	Point _p1;
	Point _p2;
	Point _p3;
	PiecewiseSpectrum _spectrum;
	Vector _orientation;
	//SurfaceProperties _props;
};

//Returns CUDA-malloc'd (managed) triangular prism whose bottom face is centered at (0,0,0)
Face* getStandardPrism();
//Returns CUDA-malloc'd (managed) square with side-lengths 2 centered at (0,0,0). Orientation is towards positive z direction
Face* getStandardPlate(PiecewiseSpectrum& const s);

class Mesh {
public:
	Mesh(Face* f, unsigned int l) : _faces(f), _numFaces(l) {};
	~Mesh() { cudaFree(_faces); };

	void transform(Transformation& const m);

	//Given ray r, returns true iff r intersects the mesh at some face. In this case, r is set to the reflected ray, s is set to the spectrum of the face, and angle is set to the angle of incidence of the intersection. 
	__device__ bool handleIntersection(Ray& r, double* s, double& angle);
	void print();
private:
	Face* _faces;
	unsigned int _numFaces;
};

class Scene {
public:
	Scene() : _meshes(NULL), _numMeshes(0), _spheres(NULL), _numSpheres(0) {};
	Scene(Mesh* m, unsigned int l) : _meshes(m), _numMeshes(l), _spheres(NULL), _numSpheres(0) {};
	Scene(Sphere* s, unsigned int l) : _meshes(NULL), _numMeshes(0), _spheres(s), _numSpheres(l) {};
	Scene(Mesh* m, unsigned int lm, Sphere* s, unsigned int ls) : _meshes(m), _numMeshes(lm), _spheres(s), _numSpheres(ls) {};
	~Scene() { cudaFree(_meshes); delete[] _spheres; };

	void transform(Transformation& const t);
	__host__ __device__ unsigned int numMeshes() { return _numMeshes; };
	__host__ __device__ Mesh* meshes() { return _meshes; };
private:
	Mesh* _meshes;
	unsigned int _numMeshes;

	Sphere* _spheres;
	unsigned int _numSpheres;
};

__global__ void traceRay(Ray* rays, Mesh* meshes, unsigned int n, double** spectrums, double* angles);