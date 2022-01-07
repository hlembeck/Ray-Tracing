#include "LightColor.cuh"
#pragma once

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
public:
	Face(Point p1, Point p2, Point p3);
	Face(Point p1, Point p2, Point p3, PiecewiseSpectrum s) : _p1(p1), _p2(p2), _p3(p3), _spectrum(s) {};
	~Face() {};

	bool onPlane(Point p);
	bool inTriangle(Point p);
	//Returns t s.t p+vt lies on face, otherwise returns -1.0
	__device__ double getIntersection(Vector v, Point p);

	__device__ bool testFace(Ray& ray, unsigned int i, PiecewiseSpectrum& s, double& angle);
protected:
	Point _p1;
	Point _p2;
	Point _p3;
	PiecewiseSpectrum _spectrum;
	Vector _orientation;
	SurfaceProperties _props;
};

class Mesh {
public:
	Mesh(Face* f, unsigned int l) : _faces(f), _numFaces(l) {};
	~Mesh() { };

	__device__ double getIntersection(Ray r);
private:
	Face* _faces;
	unsigned int _numFaces;
};

struct Scene {
	Mesh* meshes;
	unsigned int numMeshes;

	Sphere* spheres;
	unsigned int numSpheres;
};


__host__ Face* getStandardPrism();
