#include "BoundaryValueProblemSmooth.h"
#include "BoundaryValueProblemBrokenLine.h"
#include<iostream>
BoundaryValueProblemBrokenLine::BoundaryValueProblemBrokenLine() = default;

BoundaryValueProblemBrokenLine::BoundaryValueProblemBrokenLine(double kappa, complex<double> alpha, double eps) :
	BoundaryValueProblemSmooth(kappa, alpha, eps)
{
	vector<double> rho = Parameters::rho;
	vector<complex<double>> mu = Parameters::mu;
	size_t size = rho.size();
	double h = 1 / static_cast<double>(size);
	Rho = [rho, h, size](double y) {
		if (y < h / 2)
		{
			return rho[0];
		}
		if (y >= 1 - h / 2)
		{
			return rho[size - 1];
		}
		const int i = static_cast<int>((y - h / 2) / h);
		return rho[i] + (rho[i + 1] - rho[i]) / h * (y - (0.5 + i) * h);
	};
	this->Mu = [mu, h, size](double y) {
		if (y < h / 2)
		{
			return mu[0];
		}
		if (y >= 1 - h / 2)
		{
			return mu[size - 1];
		}
		const int i = static_cast<int>((y - h / 2) / h);
		return mu[i] + (mu[i + 1] - mu[i]) / h * (y - (0.5 + i) * h);
	};
}

void BoundaryValueProblemBrokenLine::CreateTheInnerSolution(size_t points, size_t size)
{
	size_t n = Parameters::rho.size();
	double h = 1.0 / n;	
	vector<complex<double>> initialValue(size);
	for (size_t i = 0; i < size; i++) initialValue[i] = complex<double>(0, 0);
	initialValue[1] = complex<double>(exp(-alpha.real()), 0);
	InnerCauchyProblemSolution.clear();
	double a = 0;
	BoundaryValueProblemSmooth::CauchyProblem(a, a+h / 2, initialValue);
	InnerCauchyProblemSolution.push_back(initialValue);
	a += h / 2;
	for (size_t i = 0; i < n - 1; i++) {
		BoundaryValueProblemSmooth::CauchyProblem(a, a + h, initialValue);
		InnerCauchyProblemSolution.push_back(initialValue);
		a += h;
	}
	BoundaryValueProblemSmooth::CauchyProblem(1 - h / 2, 1, initialValue);
	InnerCauchyProblemSolution.push_back(initialValue);
}




BoundaryValueProblemBrokenLine::~BoundaryValueProblemBrokenLine() {
}


vector<complex<double>> BoundaryValueProblemBrokenLine::CauchyProblem(double a, double b, vector<complex<double>> & InitialValue) const {
	size_t n = Parameters::rho.size();
	double h = (b - a) / n;
	BoundaryValueProblemSmooth::CauchyProblem(a, a + h / 2, InitialValue);
	a += h / 2;
	for (size_t i = 0; i < n - 1; i++) {
		BoundaryValueProblemSmooth::CauchyProblem(a, a + h, InitialValue);
		a += h;
	}
	BoundaryValueProblemSmooth::CauchyProblem(b - h / 2, b, InitialValue);
	return InitialValue;
}

