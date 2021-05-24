#pragma once
#include "Root.h"
#include<vector>
#include<complex>
#include <algorithm>
#include <cassert>

const double Pi = 3.1415926535897932384626433832795;


template <class BVP>
class SystemEvaluation
{
	vector<Root<BVP>> roots;
	double kappa;
	int counter;
	double eps;
	size_t n;
	complex<double> upperWaveField(double x) const;
public:
	SystemEvaluation();
	SystemEvaluation(double kappa, size_t n, double eps = 0.1e-5);
	SystemEvaluation(const vector<complex<double>>& Roots, double kappa,
		size_t n, double eps = 0.1e-5);
	vector<complex<double>> getRoots() const;
	int getCounter() const;
	//	~SystemEvaluation();
	vector<complex<double>>
		EvaluateTheRightPart(double a, double b, size_t N) const;
	vector<vector<complex<double>>> EvaluateTheKernel(double a, double b,
		double c, double d, size_t N);
	vector<vector<complex<double>>> EvaluateTheKernelRho(double a, double b,
		double c, double d, size_t N);
};




template<class BVP>
complex<double> SystemEvaluation<BVP>::upperWaveField(double x) const
{
	complex<double> I(0, 1);
	complex<double> result(0, 0);
	for (size_t i = 0; i < roots.size(); i++)
	{
		auto term = roots[i].u() * exp(I * roots[i].Value() * x);
		result += term;
		if (abs(term) < eps) {
			break;
		}
	}
	return I * result;
}

template<class BVP>
SystemEvaluation<BVP>::SystemEvaluation() : kappa(0), counter(0), eps(0), n(0)
{
}


template<class BVP>
SystemEvaluation<BVP>::SystemEvaluation(double kappa, size_t n, double eps) :
	kappa(kappa), counter(0), eps(eps), n(n)
{
	for (int i = 0; i < 33; i++)
	{
		auto value = sqrt(kappa * kappa / Parameters::muConst -
			Pi * Pi * (i * 2 + 1) * (1 + 2 * i) / 4);
		roots.emplace_back(kappa, value, n);
	}
	sort(roots.begin(), roots.end(), [](const Root<BVP>& a, const Root<BVP>& b)
		{
			return abs(a.Value().imag()) < abs(b.Value().imag());
		});
}

template<class BVP>
SystemEvaluation<BVP>::SystemEvaluation(const vector<complex<double>>& Roots,
	double kappa, size_t n, double eps)
	: kappa(kappa), counter(0), eps(eps), n(0)
{
	const size_t counter = Roots.size();
	for (size_t i = 0; i < counter; i++)
	{
		roots.emplace_back(kappa, Roots[i], n);//.push_back(root);
	}
	sort(roots.begin(), roots.end(), [](const Root<BVP>& a, const Root<BVP>& b)
		{
			return abs(a.Value().imag()) < abs(b.Value().imag());
		});
}

template<class BVP>
vector<complex<double>> SystemEvaluation<BVP>::getRoots() const {
	vector<complex<double>> result;
	for (size_t i = 0; i < roots.size(); i++) {
		result.push_back(roots[i].Value());
	}
	return result;
}

template <class BVP>
int SystemEvaluation<BVP>::getCounter() const
{
	return counter;
}

/*
template<class BVP>
SystemEvaluation<BVP>::~SystemEvaluation() = default;
*/
template<class BVP>
vector<complex<double>> SystemEvaluation<BVP>::EvaluateTheRightPart(double a, double b, size_t N) const
{
	assert(a >= 0);
	assert(b > a);
	complex<double> I(0, 1);
	const double h = (b - a) / N;
	auto result = vector<complex<double>>(N);
	for (size_t i = 0; i < N; i++) {
		result[i] = upperWaveField(a + h * (i + 0.5));
	}
	return result;
}


template<class BVP>
vector<vector<complex<double>>> SystemEvaluation<BVP>::EvaluateTheKernel(
	double a, double b, double c, double d, size_t N) {
	complex<double> I(0, 1);
	for (size_t i = 0; i < roots.size(); i++) roots[i].AddTheInnerSolution(N);
	const double hx = (d - c) / static_cast<double>(N);
	vector<vector<complex<double>>> Result(N);
	for (size_t i = 0; i < N; i++) {
		Result[i] = vector<complex<double>>(N);
	}
	for (size_t j = 0; j < N; j++) {
		double xi = (j + 0.5) * hx;
		for (size_t i = 0; i < N; i++) {
			complex<double> result(0, 0);
			for (size_t k = 0; k < roots.size(); k++) {
				complex<double> multiplier = exp(xi * I * roots[k].Value());
				auto a01 = roots[k].GetA01();
				auto a11 = roots[k].GetA11();
				auto a02 = roots[k].GetA02();
				auto a12 = roots[k].GetA12();
				auto b1 = roots[k].GetB1();
				auto b2 = roots[k].GetB2();
				auto r1 = I * (2.0 * a01[i] * (a11[i] - a01[i] * b2 / b1) +
					I * xi * a01[i] * a01[i]) / b1 / b1 / Parameters::mu[i] / Parameters::mu[i];
				auto r2 = I * (2.0 * a02[i] * (a12[i] - a02[i] * b2 / b1) +
					I * xi * a02[i] * a02[i]) / b1 / b1;
				result += (r1 + r2) * multiplier;
			}
			Result[j][i] = -result / static_cast<double>(N);
		}
	}
	return Result;
}

template<class BVP>
vector<vector<complex<double>>> SystemEvaluation<BVP>::EvaluateTheKernelRho(double a, double b, double c, double d, size_t N)
{
	complex<double> I(0, 1);
	for (size_t i = 0; i < roots.size(); i++) roots[i].AddTheInnerSolution(N);
	const double hx = (d - c) / static_cast<double>(N);
	vector<vector<complex<double>>> Result(N);
	for (size_t i = 0; i < N; i++) {
		Result[i] = vector<complex<double>>(N);
	}
	for (size_t j = 0; j < N; j++) {
		double xi = (j + 0.5) * hx;
		for (size_t i = 0; i < N; i++) {
			complex<double> result(0, 0);
			for (size_t k = 0; k < roots.size(); k++) {
				complex<double> multiplier = exp(xi * I * roots[k].Value());
				auto a03 = roots[k].GetA03();
				auto a13 = roots[k].GetA13();
				auto b1 = roots[k].GetB1();
				auto b2 = roots[k].GetB2();
				auto r2 = I * (2.0 * a03[i] * (a13[i] - a03[i] * b2 / b1) +
					I * xi * a03[i] * a03[i]) / b1 / b1;
				result += r2 * multiplier;
			}
			Result[j][i] = { kappa * kappa * result.real() / static_cast<double>(N),0 };
		}
	}
	return Result;
}
