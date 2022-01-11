#include "Types3D.cuh"

double piecewiseGaussian(double x, double m, double s, double t);

//Computes sum_{i=0}^n f(s_i)D_i for n = ceil((h-l)/s) - 1, s_i = l+s*i, and D_i = {s if i!= n-1 , h-s*floor((h-l)/s)) if i=n-1}
double rSum(double l, double h, double s, double* vals);

//computes Riemann sum like above for the product fg
double fg_RiemannSum(double l, double h, double s, double (*f)(double), double (*g)(double));

//returns least positive root if it exists. Otherwise returns -1.
double quadraticSolution(double a, double b, double c);

//If c = <a_0,...,a_n>, then Polynomial(n,c) represents the polynomial a_0 + a_1x + ... + a_nx^n
class Polynomial {
public:
	Polynomial(unsigned int d, double* c) : _deg(d), _coeff(c) {};
	~Polynomial() { delete[] _coeff; };

	double eval(double t);
private:
	unsigned int _deg;
	double* _coeff;
};