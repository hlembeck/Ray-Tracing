#include "Types3D.cuh"

__host__ __device__ void Vector::normalize() {
	if (_mag == 0)
		return;
	x = x / _mag;
	y = y / _mag;
	z = z / _mag;
	_mag = 1;
}

__host__ void Matrix::right_compose(Matrix& const m) {
	double* const n = m.getMatrix();
	_m[0] = _m[0] * n[0] + _m[1] * n[3] + _m[2] * n[6];
	_m[1] = _m[0] * n[1] + _m[1] * n[4] + _m[2] * n[7];
	_m[2] = _m[0] * n[2] + _m[1] * n[5] + _m[2] * n[8];

	_m[3] = _m[3] * n[0] + _m[4] * n[3] + _m[5] * n[6];
	_m[4] = _m[3] * n[1] + _m[4] * n[4] + _m[5] * n[7];
	_m[5] = _m[3] * n[2] + _m[4] * n[5] + _m[5] * n[8];

	_m[6] = _m[6] * n[0] + _m[7] * n[3] + _m[8] * n[6];
	_m[7] = _m[6] * n[1] + _m[7] * n[4] + _m[8] * n[7];
	_m[8] = _m[6] * n[2] + _m[7] * n[5] + _m[8] * n[8];
}

Transformation::Transformation() {
	cudaMallocManaged((void**)&_m, sizeof(double) * 16);
	cudaMemset(_m, 0, sizeof(double) * 15);
	_m[0] = 1.0;
	_m[5] = 1.0;
	_m[10] = 1.0;
	_m[15] = 1.0;
}

Transformation::Transformation(Point& const p) {
	cudaMallocManaged((void**)&_m, sizeof(double) * 16);
	cudaMemset(_m, 0, sizeof(double) * 15);
	_m[0] = 1.0;
	_m[3] = p.x;
	_m[5] = 1.0;
	_m[7] = p.y;
	_m[10] = 1.0;
	_m[11] = p.z;
	_m[15] = 1.0;
}

Transformation::Transformation(Vector& const v) {
	cudaMallocManaged((void**)&_m, sizeof(double) * 16);
	cudaMemset(_m, 0, sizeof(double) * 15);
	_m[0] = v.x;
	_m[5] = v.y;
	_m[10] = v.z;
	_m[15] = 1.0;
}

Transformation::Transformation(Vector& const v, double a) {
	v.normalize();
	double c = cos(a);
	double s = sin(a);
	cudaMallocManaged((void**)&_m, sizeof(double) * 16);
	cudaMemset(_m, 0, sizeof(double) * 15);
	_m[0] = v.x * v.x + (1 - v.x * v.x) * c;
	_m[1] = v.x * v.y * (1 - c) - v.z * s;
	_m[2] = v.x * v.z * (1 - c) + v.y * s;

	_m[4] = v.y * v.x * (1 - c) + v.z * s;
	_m[5] = v.y * v.y + (1 - v.y * v.y) * c;
	_m[6] = v.y * v.z * (1 - c) - v.x * s;

	_m[8] = v.z * v.x * (1 - c) - v.y * s;
	_m[9] = v.z * v.y * (1 - c) + v.x * s;
	_m[10] = v.z * v.z + (1 - v.z * v.z) * c;

	_m[15] = 1;
}