#include "Camera.cuh"

void Sensor::_computePixel(double** const spectrumStack, const double* angleStack, const unsigned int index) {
	for (int i = raysPerPixel * raysPerPixel * reflectionsPerRay - 1; i >= 0; i--) {
		//printArray(spectrumStack[i], STEPS);
		//printArray(spectrumStack[i], STEPS);
		//_pixelArray[index].add(*(spectrumStack+i*sizeof(double*)));
		_pixelArray[index].add(spectrumStack[i]);
		_pixelArray[index].scale(angleStack[i]);
	}
}

void PinholeCamera::_fillRays(Ray* rays, const unsigned int x, const unsigned int y) {
	//specified as a point (a,b), implicitly (a,b,0)
	//printf("%d %d\n", x, y);
	std::pair<double, double> pixelTopLeft(_sensor._width*(-0.5+((double)x)/_sensor._pWidth), _sensor._height*(0.5-((double)y)/_sensor._pHeight));
	double pLength = _sensor._width / _sensor._pWidth;
	for (unsigned int i = 0; i < raysPerPixel; i++) {
		for (unsigned int j = 0; j < raysPerPixel; j++) {
			rays[i * raysPerPixel + j].origin = Point(pixelTopLeft.first + pLength * ((double)j + 0.5) / raysPerPixel, pixelTopLeft.second + pLength * ((double)i - 0.5) / raysPerPixel, 0.0);
			rays[i * raysPerPixel + j].direction = Vector(_pinhole - rays[i * raysPerPixel + j].origin);
			rays[i * raysPerPixel + j].direction.normalize();
		}
	}
}

PiecewiseSpectrum* PinholeCamera::takePicture(Scene& scene) {
	const unsigned int spectrumSize = sizeof(double*) * raysPerPixel * raysPerPixel * reflectionsPerRay;
	const unsigned int angleSize = sizeof(double) * raysPerPixel * raysPerPixel *  reflectionsPerRay;

	//actually the cosine of the incidence angle.
	double *angles, **spectrums;
	Ray* rays;
	cudaMallocManaged((void**)&rays, sizeof(Ray) * raysPerPixel * raysPerPixel);
	cudaMallocManaged((void**)&spectrums, spectrumSize);
	for (unsigned int i = 0; i < raysPerPixel * raysPerPixel * reflectionsPerRay; i++) {
		cudaMallocManaged((void**)&(spectrums[i]), sizeof(double) * STEPS);
		cudaMemset(spectrums[i], 0, sizeof(double) * STEPS);
	}
	cudaMallocManaged((void**)&angles, angleSize);
	dim3 numBlocks(raysPerPixel, raysPerPixel);

	for (unsigned int i = 0; i < _sensor._pHeight; i++) {
		for (unsigned int j = 0; j < _sensor._pWidth; j++) {
			//no need to memset spectrums, rays, or angles as all will be filled by TraceRay regardless.
			_fillRays(rays, j, i);
			//printRays(rays, raysPerPixel, raysPerPixel);
			for (unsigned int k = 0; k < reflectionsPerRay; k++) {
				traceRay<<<numBlocks, 1>>>(rays, scene.meshes(), scene.numMeshes(), &spectrums[k*raysPerPixel], &angles[raysPerPixel]);
				cudaDeviceSynchronize();
			}
			//printRays(rays, raysPerPixel, raysPerPixel);
			_sensor._computePixel(spectrums, angles, i * _sensor._pWidth + j);
		}
	}
	cudaFree(rays);
	cudaFree(spectrums);
	cudaFree(angles);
	/*for (unsigned int i = 0; i < 64; i++) {
		printArray(_sensor._pixelArray[i].getPower(), STEPS);
	}*/
	return _sensor._pixelArray;
}