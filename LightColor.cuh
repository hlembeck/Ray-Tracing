#include "Types3D.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#pragma once

//in nm
constexpr unsigned short MIN_WAVELENGTH = 350;
constexpr unsigned short MAX_WAVELENGTH = 750;
constexpr unsigned char STEP_LENGTH = 10;
constexpr unsigned short STEPS = 40;

//Describes a spectrum that is piecewise constant with (MAX_WAVELENGTH-MIN_WAVELENGTH)/_step pieces. Each piece i emits power _power[i].
class PiecewiseSpectrum {
public:
	PiecewiseSpectrum() : _power(new double[STEPS]{ 0 }) {};
	PiecewiseSpectrum(double* p) : _power(p) { };
	~PiecewiseSpectrum() { delete[] _power; };

	double* getPower() const { return _power; };
	void add(PiecewiseSpectrum& const s) { for (unsigned int i = 0; i < STEPS; i++) { _power[i] = _power[i] + s.getPower()[i] }  };
	void scale(double angle) { for (unsigned int i = 0; i < STEPS; i++) { _power[i] = _power[i] * angle; } };
private:
	//note: angle is cosine of angle, not angle itself.
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
};