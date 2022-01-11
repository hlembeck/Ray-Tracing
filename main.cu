#include "Camera.cuh"

//To generate png files (https://github.com/lvandeve/lodepng)
#include "lodepng.h"

const unsigned short WIDTH = 4;
const unsigned short HEIGHT = 4;

std::vector<unsigned char> array2vector(RGBa* arr, unsigned int len) {
	std::vector<unsigned char> ret;
	for (unsigned int i = 0; i < len*4; i++) {
		ret.push_back(arr[i].r());
		ret.push_back(arr[i].g());
		ret.push_back(arr[i].b());
		ret.push_back(arr[i].a());
	}
	return ret;
}

int main() {
	cudaSetDevice(0);
	//Face* faces = getStandardPrism();
	//Mesh m1(faces, 8);
	PiecewiseSpectrum sunSpectrum = sun();
	Face* plate = getStandardPlate(sunSpectrum);
	Mesh m2(plate, 2);
	Mesh* meshes;
	cudaMallocManaged((void**)&meshes, sizeof(Mesh));
	meshes[0] = m2;
	Scene scene{ meshes,1 };
	Point p(0, 0, 3);
	Transformation t(p);
	scene.transform(t);

	Sensor sensor(WIDTH, HEIGHT, 2, 2);
	Point pinhole(0, 0, 1);
	PinholeCamera camera(sensor, pinhole);
	PiecewiseSpectrum* s = camera.takePicture(scene);
	/*for (unsigned int i = 0; i < 64; i++) {
		printArray(s[i].getPower(), STEPS);
	}*/
	RGBa* rgb = develop(s, WIDTH*HEIGHT);
	std::vector<unsigned  char> vec = array2vector(rgb, WIDTH * HEIGHT);

	unsigned error;
	if ((error= lodepng::encode("image.png", vec, WIDTH, HEIGHT)))
		printf("encoder error %s\n", lodepng_error_text(error));

	return 0;
}