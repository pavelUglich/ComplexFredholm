#include "IterativeProcess.h"



vector<complex<double>> iterative_process::_marching(const vector<complex<double>>& a) {
	vector<complex<double>> x(_size);
	vector<complex<double>> xi(_size + 1);
	vector<complex<double>> eta(_size + 1);
	xi[0] = 0;
	eta[0] = 0;
	for (size_t i = 0; i < _size; i++) {
		xi[i + 1] = _b[i] / (a[i] - _c[i] * xi[i]);
		eta[i + 1] = (_c[i] * eta[i] - _right_part[i]) / (a[i] - _c[i] * xi[i]);
	}
	x[_size - 1] = eta[_size];
	for (size_t i = 1; i < _size; i++)
		x[_size - i - 1] = xi[_size - i] * x[_size - i] + eta[_size - i];
	return x;
}

double iterative_process::_residual(vector<complex<double>> a, double alpha) {
	for (size_t i = 0; i < _size; i++)
		a[i] = -a[i] - alpha;
	_solution = _marching(a);//
	double nz = norm(_solution);
	vector<complex<double>> pz(_p1.size());
	for (size_t i = 0; i < _p2.size() - 1; i++) pz[i] = _p1[i] * _solution[i] + _p2[i] * _solution[i + 1];
	pz[_p1.size() - 1] = _p1[_p1.size() - 1] * _solution[_p1.size() - 1];
	for (size_t i = 0; i < pz.size(); i++) pz[i] -= _qtu[i];
	double npz = norm(pz);
	return step * step * npz * npz - (delta + h * nz) * (delta + h * nz);
}

void iterative_process::iterations_run() {
	//double AlphaS, AlphaN, ss, sn;
	auto alpha_s = _alpha;
	auto alpha_n = _alpha * 0.5;
	auto ss = _residual(_a, alpha_s);
	auto sn = _residual(_a, alpha_n);
	for (int i = 0; i < _iterations; i++) {
		double alpha_ = alpha_n / (1 - (1 / alpha_s) * (alpha_s - alpha_n) * sn / (sn - ss));
		ss = sn;
		sn = _residual(_a, alpha_);
		alpha_s = alpha_n;
		alpha_n = alpha_;
		if (abs(alpha_s - alpha_n) < eps) break;
	}
}

iterative_process::iterative_process(const vector<complex<double>>& p1,
	const vector<complex<double>>& p2,
	const vector<complex<double>>& RightPart,
	const vector<complex<double>>& Qtu, double Alpha, double step, double h,
	double delta, double eps, int iterations) : _p1(p1),
	_p2(p2), _qtu(Qtu), _right_part(RightPart), _alpha(Alpha), step(step), h(h),
	delta(delta), _iterations(iterations), _size(RightPart.size()), eps(eps) {
	tridiag();
	iterations_run();
}

void iterative_process::tridiag() {
	_a.resize(_size);
	_b.resize(_size);
	_c.resize(_size);
	const auto size = _p1.size();
	_a[0] = abs(_p1[0]) * abs(_p1[0]);
	//_a[_size - 1] = abs(_p1[_size - 1]) * abs(_p1[_size - 1]) + abs(_p2[_size - 1]) * abs(_p2[_size - 1]);
	if (_p2.size() == _p1.size())
	{
		_a[size - 1] = abs(_p2[size - 1]) * abs(_p2[size - 1]);
	}
	else
	{
		_a[size - 1] = abs(_p2[size - 2]) * abs(_p2[size - 2]) + abs(_p1[size - 1]) * abs(_p1[size - 1]);
	}
	_b[0] = conj(_p1[0]) * _p2[0];
	_b[size - 1] = 0;
	//#pragma omp parallel for
	for (size_t i = 1; i < size - 1; i++) {
		_a[i] = abs(_p2[i - 1]) * abs(_p2[i - 1]) + abs(_p1[i]) * abs(_p1[i]);
		_b[i] = conj(_p1[i]) * _p2[i];
	}
	_c[0] = 0;
	for (size_t i = 1; i < size; i++) _c[i] = conj(_b[i - 1]);
}


