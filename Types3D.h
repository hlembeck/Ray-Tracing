#include <math.h>
#include <stdio.h>
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

	Vector operator+ (Vector v);
	void normalize();
	Vector cross(Vector v);
	double dot(Vector v);
	double mag();
	void print();
private:
	double _mag;
};