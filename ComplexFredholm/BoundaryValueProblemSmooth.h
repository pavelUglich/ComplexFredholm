#pragma once
#include "BoundaryValueProblem.h"
#include<complex>
#include<functional>
#include"Parameters.h"


//function<complex<double>(double)> Rho;


/**
 * \brief Таблица Бутчера
 */
const vector<vector<double>> ButcherTableau = {
	{ 0.25,0.25 },
{ 0.375, 0.09375, 0.28125 },
{ 12 / 13.0, 1932 / 2197.0, -7200 / 2197.0, 7296 / 2197.0 },
{ 1, 439 / 216.0, -8, 3680 / 513.0, -845 / 4104.0 },
{ 0.5, -8 / 27.0, 2, -3544 / 2565.0, 1859 / 4104.0, -11 / 40.0 },
{ 0, 16 / 135.0, 0, 6656 / 12825.0, 28561 / 56430.0, -9 / 50.0, 2 / 55.0 },
{ 0, 25 / 216.0, 0, 1408 / 2565.0, 2197 / 4104.0, -0.2, 0 }
};

double Norma(const vector <complex<double>> & v);


/**
 * \brief Краевая задача для однородного упругого слоя 
 * (решение через приближенное решение краевой задачи)
 */
class BoundaryValueProblemSmooth :
	public BoundaryValueProblem {
	/**
	 * \brief правые части системы дифференциальных уравнений
	 */
	vector<function<complex<double>(complex<double>, double, complex<double>, const vector<complex<double>>&)>>
	EquationsSystem;
	vector<complex<double>> Evaluate(double, const vector<complex<double>> &) const;
	double EmbeddedFormula(double x, double step, vector<complex<double>> & InitialValue) const;
	//double kappa;
protected:
	function<complex<double>(double)> Mu;
	function<double(double)> Rho;
	void InitializeTheSystem();
	virtual vector<complex<double>> CauchyProblem(double, double, vector<complex<double>> &) const;

public:
	vector<complex<double>> CauchyProblemSolutions(double y, size_t size = 2) const override;
	complex<double> IntegrandDenumerator(size_t size) override;
    complex<double> waveField() const override;

	BoundaryValueProblemSmooth();
	~BoundaryValueProblemSmooth();
	BoundaryValueProblemSmooth(double kappa, complex<double> alpha, double eps = 0.1e-5/*, double nu = 0.3*/);
};

