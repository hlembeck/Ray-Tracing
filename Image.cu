#include "Image.cuh"

XYZ* Image::getXYZ() {
	XYZ* xyz = new XYZ[_len];
	for (unsigned int i = 0; i < _len; i++) {
		xyz[i] = _image[i].toXYZ();
	}
	return xyz;
}

RGBa* Image::getRGB() {
	RGBa* rgb = new RGBa[_len];
	XYZ temp = _image[0].toXYZ();
	for (unsigned int i = 1; i < _len; i++) {
		rgb[i-1] = temp.toRGB();
		temp = _image[i].toXYZ();
	}
	rgb[_len - 1] = temp.toRGB();
	return rgb;
}