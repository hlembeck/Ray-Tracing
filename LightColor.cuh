#include "Analysis.cuh"

//in nm
constexpr unsigned short MIN_WAVELENGTH = 350;
constexpr unsigned short MAX_WAVELENGTH = 750;
constexpr unsigned char STEP_LENGTH = 10;
constexpr unsigned short STEPS = 40;

//Unfortunately named, this is actually the square root of the number of rays per pixel; rays are emitted in a grid pattern from the pixel, whose dimensions are raysPerPixel^2
const int raysPerPixel = 1;

//Number of times a ray is to be reflected
const int reflectionsPerRay = 1;

void printArray(double* const arr, unsigned int const len);

class RGBa {
public:
	RGBa(unsigned char r, unsigned char g, unsigned char b, unsigned char a) : _b(b), _g(g), _r(r), _a(a) {};
	RGBa() : _b(0), _g(0), _r(0), _a(0) {};
	~RGBa() {};

	unsigned char r() { return _r; };
	unsigned char g() { return _g; };
	unsigned char b() { return _b; };
	unsigned char a() { return _a; };
private:
	unsigned char _b;
	unsigned char _g;
	unsigned char _r;
	unsigned char _a;
};

//Source: https://en.wikipedia.org/wiki/CIE_1931_color_space
double colorMatch_X(double l);
double colorMatch_Y(double l);
double colorMatch_Z(double l);

class XYZ {
public:
	XYZ() : _X(0), _Y(0), _Z(0) {};
	XYZ(double x, double y, double z) : _X(x), _Y(y), _Z(z) {};
	~XYZ() {};

	RGBa toRGB();
private:
	double _X;
	double _Y;
	double _Z;
};


//Describes a spectrum that is piecewise constant with (MAX_WAVELENGTH-MIN_WAVELENGTH)/_step pieces. Each piece i emits power _power[i].
class PiecewiseSpectrum {
public:
	PiecewiseSpectrum() : _power(NULL) {
		cudaMallocManaged((void**)&_power, sizeof(double) * STEPS);
		cudaMemset(_power, 0, sizeof(double) * STEPS);
	};
	PiecewiseSpectrum(double* p) : _power(p) { };
	~PiecewiseSpectrum() { cudaFree(_power); };

	__host__ __device__ double* getPower() const { return _power; };
	void add(double* const s) {
		for (unsigned int i = 0; i < STEPS; i++) {
			//printf("%f", _power[i]);
			_power[i] = _power[i] + s[i];
		}
	};
	//note: angle is cosine of angle, not angle itself.
	void scale(double angle) { for (unsigned int i = 0; i < STEPS; i++) { _power[i] = _power[i] * angle; } };
	XYZ toXYZ();
private:
	double* _power;
};

const PiecewiseSpectrum sun();

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

struct Ray {
	Point origin;
	Vector direction;
};

void printRays(Ray* rays, unsigned int w, unsigned int h);

XYZ* getXYZ(PiecewiseSpectrum* s, unsigned int l);
RGBa* getRGB(XYZ* xyz, unsigned int l);

RGBa* develop(PiecewiseSpectrum* s, unsigned int l);