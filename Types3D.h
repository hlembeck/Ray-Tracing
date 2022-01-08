#include <math.h>
#include <stdio.h>
#include <utility>
#pragma once

class Vector;

class Point {
public:
	Point() : x(0.0), y(0.0), z(0.0) {};
	Point(double a, double b, double c) : x(a), y(b), z(c) {};
	~Point() {};

	Point operator + (Point& p);
	Point operator - (Point& p);
	bool equals(Point p);
	Point translate(Vector v, double t);
	double x;
	double y;
	double z;
};

class Vector {
public:
	Vector(double a, double b, double c) : x(a), y(b), z(c), _mag(sqrt(a* a + b * b + c * c)) {};
	Vector(Point p) : x(p.x), y(p.y), z(p.z), _mag(sqrt(x* x + y * y + z * z)) {};
	Vector() : x(0.0), y(0.0), z(0.0), _mag(0.0) {};
	~Vector() {};
	double x;
	double y;
	double z;

	Vector operator+ (const Vector& v);
	Vector operator- (const Vector& v);
	Vector scale(const double t);
	void normalize();
	Vector cross(Vector v);
	double dot(Vector v);
	double mag();
	void print();
private:
	double _mag;
};


class Matrix {
public:
	Matrix(double* mat) : _m(mat) {};
	~Matrix() { delete[] _m };

	Vector left_transform(Vector& const v) { return Vector(_m[0] * v.x + _m[1] * v.y + _m[2] * v.z, _m[3] * v.x + _m[4] * v.y + _m[5] * v.z, _m[6] * v.x + _m[7] * v.y + _m[8] * v.z); };
private:
	double* _m;
};