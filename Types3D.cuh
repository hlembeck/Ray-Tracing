#include <math.h>
#include <stdio.h>
#include <string.h>
#include <utility>
#include <stdlib.h>
#include <limits>
#include <vector>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

class Vector;

class Point {
public:
	__host__ __device__ Point() : x(0.0), y(0.0), z(0.0) {};
	__host__ __device__ Point(double a, double b, double c) : x(a), y(b), z(c) {};
	__host__ __device__ Point(Vector& const v);
	__host__ __device__ ~Point() {};

	__host__ __device__ Point operator + (Point& p) { return Point(x + p.x, y + p.y, z + p.z); };
	__host__ __device__ Point operator - (Point& p) { return Point(x - p.x, y - p.y, z - p.z); };
	__host__ __device__ bool equals(Point p) { return (x == p.x ? (y == p.y ? (z == p.z ? true : false) : false) : false); };
	__host__ __device__ Point translate(Vector v, double t);
	__host__ __device__ void print() { printf("[%f %f %f]", x, y, z); };
	double x;
	double y;
	double z;
};

class Vector {
public:
	__host__ __device__ Vector(double a, double b, double c) : x(a), y(b), z(c), _mag(sqrt(a* a + b * b + c * c)) {};
	__host__ __device__ Vector(Point p) : x(p.x), y(p.y), z(p.z), _mag(sqrt(x* x + y * y + z * z)) {};
	__host__ __device__ Vector() : x(0.0), y(0.0), z(0.0), _mag(0.0) {};
	__host__ __device__ ~Vector() {};
	double x;
	double y;
	double z;

	__host__ __device__ Vector operator+ (const Vector& v) { return Vector(x + v.x, y + v.y, z + v.z); };
	__host__ __device__ Vector operator- (const Vector& v) { return Vector(x - v.x, y - v.y, z - v.z); };
	__host__ __device__ Vector scale(const double t) { return Vector(x * t, y * t, z * t); };
	__host__ __device__ void normalize();
	__host__ __device__ Vector cross(Vector v) { return Vector(y * v.z - z * v.y, v.x * z - x * v.z, x * v.y - y * v.x); };
	__host__ __device__ double dot(Vector v) { return x * v.x + y * v.y + z * v.z; };
	__host__ __device__ double mag() { return _mag; };
	void print() { printf("[%f %f %f]", x, y, z); };
private:
	double _mag;
};

inline Point::Point(Vector& const v) : x(v.x), y(v.y), z(v.z) {}
inline Point Point::translate(Vector v, double t) { return Point(x + v.x * t, y + v.y * t, z + v.z * t); }

class Matrix {
public:
	Matrix(double* mat) : _m(mat) {};
	~Matrix() { delete[] _m; };

	double* getMatrix() { return _m; };

	//Given vector v, computes Mv (vector is left-transformed).
	Vector left_transform(Vector& const v) { return Vector(_m[0] * v.x + _m[1] * v.y + _m[2] * v.z, _m[3] * v.x + _m[4] * v.y + _m[5] * v.z, _m[6] * v.x + _m[7] * v.y + _m[8] * v.z); };
	//Given matrix N, computes MN (M is right-composed with N).
	void right_compose(Matrix& const m);
private:
	double* _m;
};

class Transformation {
public:
	//Identity transformation
	Transformation();
	Transformation(double* m) : _m(m) {};
	//Defines the transformation that translates points with (0,0,0)->p
	Transformation(Point& const p);
	//Defines the transformation given by (x,y,z)->(v.x*x,v.y*y,v.z*z)
	Transformation(Vector& const v);
	//Defines the transformation that rotates some angle about the axis given by v.
	Transformation(Vector& const v, double a);
	~Transformation() { cudaFree(_m); };

	double* getTransformation() { return _m; };

	void print() { printf("test"); printf("%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n", _m[0], _m[1], _m[2], _m[3], _m[4], _m[5], _m[6], _m[7], _m[8], _m[9], _m[10], _m[11], _m[12], _m[13], _m[14], _m[15]); };

	void transform(Point& p) {
		p.x = _m[0] * p.x + _m[1] * p.y + _m[2] * p.z + _m[3];
		p.y = _m[4] * p.x + _m[5] * p.y + _m[6] * p.z + _m[7];
		p.z = _m[8] * p.x + _m[9] * p.y + _m[10] * p.z + _m[11];
	};
private:
	double* _m;
};