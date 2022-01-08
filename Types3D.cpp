#include "Types3D.h"

Point Point::operator +(Point& p) {
	return Point(x + p.x, y + p.y, z + p.z);
}

Point Point::operator -(Point& p) {
	return Point(x - p.x, y - p.y, z - p.z);
}

bool Point::equals(Point p) {
	return (x == p.x ? (y == p.y ? (z == p.z ? true : false) : false) : false);
}

Point Point::translate(Vector v, double t) {
	return Point(x + v.x * t, y + v.y * t, z + v.z * t);
}

Vector Vector::cross(Vector v) {
	return Vector(y * v.z - z * v.y, v.x * z - x * v.z, x * v.y - y * v.x);
}

double Vector::dot(Vector v) {
	return x * v.x + y * v.y + z * v.z;
}

void Vector::normalize() {
	if (_mag == 0)
		return;
	x = x / _mag;
	y = y / _mag;
	z = z / _mag;
	_mag = 1;
}

double Vector::mag() {
	return _mag;
}

Vector Vector::operator+(const Vector& v) {
	return Vector(x + v.x, y + v.y, z + v.z);
}

Vector Vector::operator-(const Vector& v) {
	return Vector(x - v.x, y - v.y, z - v.z);
}

Vector Vector::scale(const double t) {
	return Vector(x * t, y * t, z * t);
}

void Vector::print() {
	printf("[%f %f %f]\n", x, y, z);
}