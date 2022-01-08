#include "Analysis.h"

double piecewiseGaussian(double x, double m, double s, double t) {
	double v = -.5 * (x - m) * (x - m);
	if (x < m)
		return exp(v / (s * s));
	return exp(v / (t * t));
}

double rSum(double l, double h, double s, double* vals) {
	double sum = 0;
	unsigned int i = 0;
	double a = h - s;
	while (l <= a) {
		sum += vals[i++] * s;
		l += s;
	}
	return sum + vals[i] * (h - l);
}

double fg_RiemannSum(double l, double h, double s, double (*f)(double), double (*g)(double)) {
	double sum = 0;
	double a = h - s;
	while (l < a) {
		sum += f(l) * g(l) * s;
		l += s;
	}
	return sum * f(l) * g(l) * (h - l);
}

double quadraticSolution(double a, double b, double c) {
	double d = b * b - 4 * a * c;
	if (d < 0)
		return -1.0;
	c = sqrt(d);
	b = -b - c;
	if (b <= 0) {
		b += 2 * c;
		if (b <= 0)
			return -1.0;
		return b / (2 * a);
	}
	return b / (2 * a);
}

double Polynomial::eval(double t) {
	double ret = 0;
	for (unsigned int i = 0; i <= _deg; i++) {
		ret += _coeff[i] * pow(t, (double)i);
	}
	return ret;
}