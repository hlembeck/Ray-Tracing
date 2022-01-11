#include "LightColor.cuh"

class Image {
public:
	Image(PiecewiseSpectrum* s, unsigned int l) : _image(s), _len(l) { };
	~Image() {delete[] _image; };

	XYZ* getXYZ();
	RGBa* getRGB();
private:
	PiecewiseSpectrum* _image;
	const unsigned int _len;
};