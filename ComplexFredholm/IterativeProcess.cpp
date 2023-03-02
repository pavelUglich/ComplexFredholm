#include "IterativeProcess.h"



vector<complex<double>> IterativeProcess::Marching(const vector<complex<double>> & a) {
	vector<complex<double>> x(size);
	vector<complex<double>> xi(size + 1);
	vector<complex<double>> eta(size + 1);
	xi[0] = 0;
	eta[0] = 0;
	for (size_t i = 0; i < size; i++) {
		xi[i + 1] = b[i] / (a[i] - c[i] * xi[i]);
		eta[i + 1] = (c[i] * eta[i] - RightPart[i]) / (a[i] - c[i] * xi[i]);
	}
	x[size - 1] = eta[size];
	for (size_t i = 1; i < size; i++)
		x[size - i - 1] = xi[size - i] * x[size - i] + eta[size - i];
	return x;
}

double IterativeProcess::Residual(vector<complex<double>> a, double alpha) {
	for (size_t i = 0; i < size; i++)
		a[i] = -a[i] - alpha;
	Solution = Marching(a);//
	double nz = norm(Solution);
	vector<complex<double>> pz(size);
	for (size_t i = 0; i < size - 1; i++) pz[i] = p1[i] * Solution[i] + p2[i] * Solution[i + 1];
	pz[size - 1] = p1[size - 1] * Solution[size - 1];
	for (size_t i = 0; i < size; i++) pz[i] -= Qtu[i];
	double npz = norm(pz);
	return step*step*npz*npz - (delta + h*nz)*(delta + h*nz);
}

void IterativeProcess::IterationsRun() {
	double AlphaS, AlphaN, ss, sn;
	AlphaS = alpha;
	AlphaN = alpha*0.5;
	ss = Residual(a, AlphaS);
	sn = Residual(a, AlphaN);
	for (int i = 0; i < Iterations; i++) {
		double alpha_ = AlphaN / (1 - (1 / AlphaS)*(AlphaS - AlphaN)*sn / (sn - ss));
		ss = sn;
		sn = Residual(a, alpha_);
		AlphaS = AlphaN;
		AlphaN = alpha_;
		if (abs(AlphaS - AlphaN) < eps) break;
	}
}

void IterativeProcess::tridiag() {
	a.resize(size);
	b.resize(size);
	c.resize(size);
	a[0] = abs(p1[0]) * abs(p1[0]);
	a[size - 1] = abs(p1[size - 1]) * abs(p1[size - 1]) + abs(p2[size - 1]) *abs(p2[size - 1]);
	b[0] = conj(p1[0]) * p2[0];
	b[size - 1] = 0;
	//#pragma omp parallel for
	for (size_t i = 1; i < size - 1; i++) {
		a[i] = abs(p2[i - 1]) * abs(p2[i - 1]) + abs(p1[i]) * abs(p1[i]);
		b[i] = conj(p1[i]) * p2[i];
	}
	c[0] = 0;
	for (size_t i = 1; i < size; i++) c[i] = conj(b[i - 1]);
}


