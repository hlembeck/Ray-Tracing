#include "Camera.cuh"

void Sensor::_computePixel(PiecewiseSpectrum* const spectrumStack, const double* angleStack, const unsigned int index) {
	for (unsigned int i = raysPerPixel * raysPerPixel - 1; i >= 0; i--) {
		_pixelArray[index].add(spectrumStack[i]);
		_pixelArray[index].scale(angleStack[i]);
	}
}

__host__ void PinholeCamera::_fillRays(Ray* rays, const unsigned int x, const unsigned int y) {
	//specified as a point (x,y), implicitly (x,y,0)
	std::pair<double, double> pixelTopLeft(_sensor._width*(-0.5+((double)x)/_sensor._pWidth), _sensor._height*(0.5-((double)y)/_sensor._pHeight));
	double pLength = _sensor._width / _sensor._pWidth;
	for (unsigned int i = 0; i < raysPerPixel; i++) {
		rays[i].origin = Point(pixelTopLeft.first+pLength*((double)i+0.5)/raysPerPixel, pixelTopLeft.second+pLength*((double)i+.5)/raysPerPixel, 0.0);
		rays[i].direction = Vector(_pinhole-rays[i].origin);
		rays[i].direction.normalize();
	}
}

__global__ void TraceRay(Ray* rays, Scene& scene, PiecewiseSpectrum* spectrums, double* angles) {
	int x = blockIdx.x;
	int y = blockIdx.y;
	if (x < raysPerPixel && y<raysPerPixel) {
		for (unsigned int i = 0; i < scene.numMeshes; i++) {
			if (scene.meshes[i].handleIntersection(rays[x * raysPerPixel + y], spectrums[x * raysPerPixel + y], angles[x * raysPerPixel + y]))
				goto fin;
		}
	}
	spectrums[x * raysPerPixel + y] = PiecewiseSpectrum();
	angles[x * raysPerPixel + y] = 1.0;
	fin:
}

PiecewiseSpectrum* PinholeCamera::takePicture(double exposure, Scene& scene) {
	const unsigned int spectrumSize = sizeof(PiecewiseSpectrum) * raysPerPixel * raysPerPixel * reflectionsPerRay;
	const unsigned int angleSize = sizeof(double) * raysPerPixel * raysPerPixel *  reflectionsPerRay;

	PiecewiseSpectrum* spectrums;
	//actually the cosine of the incidence angle.
	double* angles;
	Ray* rays;
	cudaMalloc((void**)&rays, sizeof(Ray) * raysPerPixel * raysPerPixel);
	cudaMalloc((void**)&spectrums, spectrumSize);
	cudaMalloc((void**)&angles, angleSize);

	for (unsigned int i = 0; i < _sensor._pHeight; i++) {
		for (unsigned int j = 0; j < _sensor._pWidth; i++) {
			//no need to memset spectrums, rays, or angles as all will be filled by TraceRay regardless.
			_fillRays(rays, j, i);
			for (unsigned int k = 0; k < reflectionsPerRay; k++) {
				TraceRay << <1, dim3(raysPerPixel, raysPerPixel) >> > (rays, scene, &spectrums[k*raysPerPixel], &angles[raysPerPixel]);
			}


		}
	}
}