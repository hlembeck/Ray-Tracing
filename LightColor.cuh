#include "Types3D.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#pragma once

//in nm
constexpr unsigned short MIN_WAVELENGTH = 350;
constexpr unsigned short MAX_WAVELENGTH = 750;

//Describes a spectrum that is piecewise constant with (MAX_WAVELENGTH-MIN_WAVELENGTH)/_step pieces. Each piece i emits power _power[i].
class PiecewiseSpectrum {
public:
	PiecewiseSpectrum() : _step(MAX_WAVELENGTH - MIN_WAVELENGTH), _power(new double[1]{ 0 }) {};
	PiecewiseSpectrum(unsigned char s, double* p) : _step(s), _power(p) { (MAX_WAVELENGTH - MIN_WAVELENGTH) % s ? throw : 0; };
	~PiecewiseSpectrum() { delete[] _power; };

private:
	unsigned short _step;
	double* _power;
};

//Describes a "ray" of some spectrum that emits from _origin in _direction for _duration seconds.
class RayBeam {
public:

private:
	PiecewiseSpectrum _spectrum;
	Point _origin;
	Vector _direction;
	double _duration;
};

//Refractive index approximated using Cauchy's equation. The refractive index for a wavelength x (in the visual spectrum, in nm) is about r1+r2/x^2.
struct RefractiveIndex {
	double r1;
	double r2;
};

__device__ struct Ray {
	Point origin;
	Vector direction;
	PiecewiseSpectrum spectrum;
};