#include "Surface.cuh"

const int raysPerPixel = 1;

//Sensor is to exist on the plane z=0, centered at (0,0,0)
__host__ struct Sensor {
	PiecewiseSpectrum* pixelArray;
	unsigned int len;
	double width;
	double height;
};

__host__ class PinholeCamera {
public:
	PinholeCamera(Sensor s, Point p) : _sensor(s), _pinhole(p) {};
	~PinholeCamera() {};

	PiecewiseSpectrum* takePicture(double exposure, Scene s);
private:
	Sensor _sensor;
	Point _pinhole;
};