#include "BoundaryValueProblemSmooth.h"
#include"Parameters.h"
#include <cassert>
complex<double> I(0, 1);
const double Pi = 3.1415926535897932384626433832795;

/**
 * \brief Расчет правой части
 * \param t поперечная координата
 * \param vec вектор неизвестных
 * \return значения правой части
 */
vector<complex<double>> BoundaryValueProblemSmooth::Evaluate(double t, 
	const vector<complex<double>> & vec) const {
	assert(vec.size() % 2 == 0);
	vector<complex<double>> Result(vec.size());
	for (size_t i = 0; i < vec.size(); i++)
	{
		Result[i] = EquationsSystem[i](kappa, t, alpha, vec);
	}
	return Result;
}

/**
 * \brief Расчёт по вложенной формуле Рунге
 * \param x поперечная координата
 * \param step длина отрезка 
 * \param InitialValue значения неизвестных функций на левом конце отрезка
 * \return значения неизвестных функций на правом конце отрезка
 */
double BoundaryValueProblemSmooth::EmbeddedFormula(double x, double step, 
	vector<complex<double>> & InitialValue) const {
	auto y = Evaluate(x, InitialValue);
	vector<vector<complex<double>>> k;
	vector<complex<double>> kk(y.size());
	for (size_t i = 0; i < y.size(); i++) kk[i] = step * y[i];
	k.push_back(kk);
	for (int i = 0; i < 5; i++) {
		auto yh = InitialValue;
		for (int j = 0; j <= i; j++)
			for (size_t l = 0; l < yh.size(); l++)
				yh[l] += ButcherTableau[i][j + 1] * k[j][l];
		auto y2 = Evaluate(x + ButcherTableau[i][0] * step, yh);
		for (size_t i = 0; i < y.size(); i++) kk[i] = step * y2[i];
		k.push_back(kk);
	}
	vector<complex<double>> y1(y.size());
	vector<complex<double>> y2(y.size());
	for (size_t i = 0; i < y.size(); i++) {
		y1[i] = 0;
		y2[i] = 0;
		for (int j = 0; j <= 5; j++) {
			y1[i] += ButcherTableau[5][j + 1] * k[j][i];
			y2[i] += ButcherTableau[6][j + 1] * k[j][i];
		}
	}
	for (size_t i = 0; i < y1.size(); i++) y2[i] -= y1[i];
	InitialValue = y1;
	return Norma(y2) / Norma(y1);
}

/**
 * \brief Численное решение задачи Коши
 * \param a левый конец отрезка
 * \param b правый конец отрезка
 * \param InitialValue начальные условия
 * \return решение на правом конце отрезка
 */
vector<complex<double>> BoundaryValueProblemSmooth::CauchyProblem(double a, 
	double b, vector<complex<double>> & InitialValue) const {
	double x = a;
	double h = b - a;
	while (abs(b - x) > eps) {
		auto yp = InitialValue;
		double n = EmbeddedFormula(x, h, yp);
		while (n > eps) {
			yp = InitialValue;
			h = 0.5*h;
			n = EmbeddedFormula(x, h, yp);
		}
		x += h;
		for (size_t i = 0; i < yp.size(); i++) InitialValue[i] += yp[i];
	}
	return InitialValue;
}

vector<complex<double>> BoundaryValueProblemSmooth::CauchyProblemSolutions(double y, size_t size) const
{
	vector<complex<double>> initialValue(size, complex<double>(0, 0));
	initialValue[1] = complex<double>(exp(-alpha.real()), 0);
	return CauchyProblem(0, y, initialValue);
}


double /*BoundaryValueProblemSmooth::*/Norma(const vector<complex<double>> & v) {
	double sum = 0;
	for (auto i : v) sum += i.real()* i.real() + i.imag()* i.imag();
	return sqrt(sum);
}

void BoundaryValueProblemSmooth::InitializeTheSystem() {
	EquationsSystem = {
		[this](complex<double> kappa, double t, complex<double> alpha, const vector<complex<double>> & vec) { return vec[1] / Mu(t); },
		[this](complex<double> kappa, double t, complex<double> alpha, const vector<complex<double>> & vec) { return vec[0] *(alpha*alpha*  Mu(t) - kappa * kappa); },
		[this](complex<double> kappa, double t, complex<double> alpha, const vector<complex<double>> & vec) { return vec[3] / Mu(t); },
		[this](complex<double> kappa, double t, complex<double> alpha, const vector<complex<double>> & vec) { return vec[0] * Mu(t)*2.0*alpha + vec[2] * (alpha*alpha* Mu(t) - kappa * kappa); },
		[this](complex<double> kappa, double t, complex<double> alpha, const vector<complex<double>> & vec) { return vec[5] / Mu(t); },
		[this](complex<double> kappa, double t, complex<double> alpha, const vector<complex<double>> & vec) { return vec[0] * Mu(t)*2.0 + vec[2] * Mu(t)*4.0*alpha + vec[4] * (alpha*alpha* Mu(t) - kappa * kappa); }
	};
}


BoundaryValueProblemSmooth::BoundaryValueProblemSmooth()
{
}


BoundaryValueProblemSmooth::~BoundaryValueProblemSmooth()
{
}

BoundaryValueProblemSmooth::BoundaryValueProblemSmooth(double kappa, complex<double> alpha, double eps) :BoundaryValueProblem(kappa, alpha, eps) {
	this->Mu = Parameters::Mu;
	this->Rho = Parameters::Rho;
	InitializeTheSystem();
}


/**
 * \brief Знаменатель подынтегрального выражения
 * \param size размер системы уравнений
 * \return Знаменатель подынтегрального выражения
 */
complex<double> BoundaryValueProblemSmooth::IntegrandDenumerator(size_t size) {
	vector<complex<double>> y(size, complex<double>(0, 0));
	y[1] =  complex<double>(exp(-alpha.real()),0);
	CauchyProblemSolution = CauchyProblem(0, 1, y);
	if (size == 2) return CauchyProblemSolution[1];
	if (size == 4) return CauchyProblemSolution[1] / CauchyProblemSolution[3];
	if (size == 6) return CauchyProblemSolution[5];
	return 0;
}

complex<double> BoundaryValueProblemSmooth::waveField() const
{
	vector<complex<double>> y(2);
	y[0] = { 0, 0 };
	y[1] = { exp(-alpha.real()), 0 };
	auto cauchyProblemSolution = CauchyProblem(0, 1, y);
	return cauchyProblemSolution[0] / cauchyProblemSolution[1];
}

