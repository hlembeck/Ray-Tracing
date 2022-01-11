#include "Surface.cuh"

//Sensor is to exist on the plane z=0, centered at (0,0,0)
class Sensor {
	friend class PinholeCamera;

public:
	Sensor(const unsigned int pw, const unsigned int ph, const double w, const double h) : _pWidth(pw), _pHeight(ph), _width(w), _height(h), _pixelArray(NULL) {};
	~Sensor() { delete[] _pixelArray; };
private:
	void _computePixel(double** const spectrumStack, const double* angleStack, const unsigned int index);
	PiecewiseSpectrum* _pixelArray;
	//Assumed that the ratio _pWidth/_pHeight is equal to _width/_height, so that pixels are square
	const unsigned int _pWidth;
	const unsigned int _pHeight;
	const double _width;
	const double _height;
};

//struct Sensor {
//	PiecewiseSpectrum* _pixelArray;
//	//Assumed that the ratio _pWidth/_pHeight is equal to _width/_height, so that pixels are square
//	const unsigned int _pWidth;
//	const unsigned int _pHeight;
//	const double _width;
//	const double _height;
//};


class PinholeCamera {
public:
	PinholeCamera(Sensor s, Point& p) : _sensor(s), _pinhole(p) { _sensor._pixelArray = new PiecewiseSpectrum[_sensor._pWidth * _sensor._pHeight]; };
	~PinholeCamera() {};

	__host__ void _fillRays(Ray* rays, const unsigned int x, const unsigned int y);
	PiecewiseSpectrum* takePicture(Scene& scene);
private:
	Sensor _sensor;
	Point _pinhole;
};