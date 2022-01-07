#include "Camera.cuh"

__global__ void TraceRay(Ray* rays, unsigned int l,  Scene scene) {
	int x = blockIdx.x;
	double t;
	double sign = 0;
	if (x < l) {
		for (unsigned int i = 0; i < scene.numMeshes; i++) {
			if ((t=scene.meshes[i].getIntersection(rays[x])) > 0) {
				rays[x].origin = rays[x].origin.translate(rays[x].direction, t);
				break;
			}
			if ((t = scene.meshes[i].getIntersection(rays[x])) < 0) {
				rays[x].origin = rays[x].origin.translate(rays[x].direction, t);
				break;
			}
		}
	}
}

PiecewiseSpectrum* PinholeCamera::takePicture(double exposure, Scene s) {
	for (unsigned int i = 0; i < _sensor.len; i++) {
		Ray* rays;
		cudaMalloc((void**)&rays, sizeof(Ray) * raysPerPixel);
		TraceRay <<<1, raysPerPixel>>>(rays, s);
	}
}