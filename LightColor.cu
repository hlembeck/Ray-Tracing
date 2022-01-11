#include "LightColor.cuh"

void printArray(double* const arr, unsigned int const len) {
	for (unsigned int i = 0; i < len - 1; i++) {
		printf("%f, ", arr[i]);
	}
	printf("%f\n", arr[len - 1]);
}

//x: wavelength (nm)
double redSpectrum(double x) {
	return .25 * exp(.5 * (x - 620) * (x - 620) / 900.0);
}

double greenSpectrum(double x) {
	return .2 * exp(.5 * (x - 535.0) * (x - 535.0) / 900.0);
}

double blueSpectrum(double x) {
	return .45 * exp(-.5 * (x - 460.0) * (x - 460.0) / 900.0);
}

double colorMatch_X(double l) {
	return 1.056 * piecewiseGaussian(l, 599.8, 37.9, 31.0) + .362 * piecewiseGaussian(l, 442.0, 16.0, 26.7) - .065 * piecewiseGaussian(l, 501.1, 20.4, 26.2);
}

double colorMatch_Y(double l) {
	return 0.821 * piecewiseGaussian(l, 568.8, 46.9, 40.5) + .286 * piecewiseGaussian(l, 530.9, 16.3, 31.1);
}

double colorMatch_Z(double l) {
	return 1.217 * piecewiseGaussian(l, 437.0, 11.8, 36.0) + .681 * piecewiseGaussian(l, 459.0, 26.0, 13.8);
}

RGBa XYZ::toRGB() {
	double* mat = new double[9];

	mat[0] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &redSpectrum, &colorMatch_X);
	mat[1] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &redSpectrum, &colorMatch_Y);
	mat[2] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &redSpectrum, &colorMatch_Z);

	mat[3] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &greenSpectrum, &colorMatch_X);
	mat[4] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &greenSpectrum, &colorMatch_Y);
	mat[5] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &greenSpectrum, &colorMatch_Z);

	mat[6] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &blueSpectrum, &colorMatch_X);
	mat[7] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &blueSpectrum, &colorMatch_Y);
	mat[8] = fg_RiemannSum(MIN_WAVELENGTH, MAX_WAVELENGTH, 1, &blueSpectrum, &colorMatch_Z);

	Matrix m(mat);
	Vector v(_X, _Y, _Z);
	Vector rgb = m.left_transform(v);
	return RGBa(rgb.x, rgb.y, rgb.z, 255);
}

XYZ PiecewiseSpectrum::toXYZ() {
	double* tVals = new double[STEPS];
	double x, y, z;
	for (unsigned short i = 0; i < STEPS; i++) {
		tVals[i] = _power[i] * colorMatch_X(MIN_WAVELENGTH + STEP_LENGTH * i);
	}
	x = rSum(MIN_WAVELENGTH, MAX_WAVELENGTH, STEP_LENGTH, tVals);

	for (unsigned short i = 0; i < STEPS; i++) {
		tVals[i] = _power[i] * colorMatch_Y(MIN_WAVELENGTH + STEP_LENGTH * i);
	}
	y = rSum(MIN_WAVELENGTH, MAX_WAVELENGTH, STEP_LENGTH, tVals);

	for (unsigned short i = 0; i < STEPS; i++) {
		tVals[i] = _power[i] * colorMatch_Z(MIN_WAVELENGTH + STEP_LENGTH * i);
	}
	z = rSum(MIN_WAVELENGTH, MAX_WAVELENGTH, STEP_LENGTH, tVals);

	return XYZ(x, y, z);
}

void printRays(Ray* rays, unsigned int w, unsigned int h) {
	printf("Origins:\n");
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; j++) {
			rays[i * w + j].origin.print();
			printf(" ");
		}
		printf("\n");
	}
	printf("Directions:\n");
	for (unsigned int i = 0; i < h; i++) {
		for (unsigned int j = 0; j < w; j++) {
			rays[i * w + j].direction.print();
			printf(" ");
		}
		printf("\n");
	}
}

XYZ* getXYZ(PiecewiseSpectrum* s, unsigned int l) {
	XYZ* xyz = new XYZ[l];
	for (unsigned int i = 0; i < l; i++) {
		xyz[i] = s[i].toXYZ();
	}
	return xyz;
}

RGBa* getRGB(XYZ* xyz, unsigned int l) {
	RGBa* rgb = new RGBa[l];
	for (unsigned int i = 0; i < l; i++) {
		rgb[i] = xyz[i].toRGB();
	}
	return rgb;
}

RGBa* develop(PiecewiseSpectrum* s, unsigned int l) {
	if (l == 0)
		return NULL;
	RGBa* rgb = new RGBa[l];
	XYZ* xyz = new XYZ[l];
	xyz[0] = s[0].toXYZ();
	for (unsigned int i = 1; i < l; i++) {
		xyz[i] = s[i].toXYZ();
		rgb[i - 1] = xyz[i - 1].toRGB();
	}
	rgb[l - 1] = xyz[l - 1].toRGB();
	return rgb;
}

const PiecewiseSpectrum sun() {
	double* vals;
	cudaMallocManaged((void**)&vals, sizeof(double) * STEPS);
	int peak = STEPS / 3;
	for (int i = 0; i < peak; i++) {
		vals[i] = 1.5 + 0.25 * (i / (double)peak);
	}
	for (int i = 0; i < STEPS-peak; i++) {
		vals[i+peak] = 1.75 - 0.35 * (i / (double)(STEPS-peak));
	}
	return PiecewiseSpectrum(vals);
}