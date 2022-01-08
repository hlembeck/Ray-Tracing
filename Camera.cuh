#include "Surface.cuh"

//Unfortunately named, this is actually the square root of the number of rays per pixel; rays are emitted in a grid pattern from the pixel, whose dimensions are raysPerPixel^2
const int raysPerPixel = 2;

//Number of times a ray is to be reflected
const int reflectionsPerRay = 5;

//Sensor is to exist on the plane z=0, centered at (0,0,0)
__host__ class Sensor {
	friend class PinholeCamera;
private:
	void _computePixel(PiecewiseSpectrum* const spectrumStack, const double* angleStack, const unsigned int index);
	PiecewiseSpectrum* _pixelArray;
	//Assumed that the ratio _pWidth/_pHeight is equal to _width/_height, so that pixels are square
	unsigned int _pWidth;
	unsigned int _pHeight;
	double _width;
	double _height;
};


__host__ class PinholeCamera {
public:
	PinholeCamera(Sensor s, Point p) : _sensor(s), _pinhole(p) {};
	~PinholeCamera() {};

	__host__ void _fillRays(Ray* rays, const unsigned int x, const unsigned int y);
	PiecewiseSpectrum* takePicture(Scene& scene);
private:
	Sensor _sensor;
	Point _pinhole;
};