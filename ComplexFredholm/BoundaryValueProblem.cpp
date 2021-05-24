#include "BoundaryValueProblem.h"
#include"Parameters.h"
#include<vector>

using namespace std;

complex<double> BoundaryValueProblem::IntegrandDenumerator(size_t size) {
	CauchyProblemSolution = CauchyProblemSolutions(1, size);
	if (size == 2) return CauchyProblemSolution[1];
	if (size == 4) return  CauchyProblemSolution[1] / CauchyProblemSolution[3];
	return 0.0;
}

complex<double> BoundaryValueProblem::waveField() const
{
	auto cauchyProblemSolution = CauchyProblemSolutions(1, 2);
	return cauchyProblemSolution[0] / cauchyProblemSolution[1];
}

complex<double> BoundaryValueProblem::perturb(double y) const
{
	const auto kappa1 = kappa / sqrt(Parameters::muConst);
	const auto gamma = sqrt(alpha * alpha - kappa1 * kappa1);
	auto sm = sinh(gamma * y);
	sm *= sm * alpha * alpha;
	auto cm = cosh(gamma * y);
	cm *= cm * gamma * gamma;
	auto numerator = (cm + sm) / gamma / gamma;
	cm = cosh(gamma);
	numerator /= cm * cm;
	return numerator;
}


void BoundaryValueProblem::SetAlpha(complex<double> _alpha) {
	alpha = _alpha;
}

void BoundaryValueProblem::CreateTheInnerSolution(size_t points, size_t size) {
	const double h = 1 / static_cast<double>(points);
	InnerCauchyProblemSolution.clear();
	for (size_t i = 0; i < points; i++) {
		InnerCauchyProblemSolution.push_back(CauchyProblemSolutions((i + 0.5) * h, size));
	}
}

double BoundaryValueProblem::getKappa() const
{
	return kappa;
}


///<summary>
///Метод Ньютона для отыскания комплексных корней дисперсионного уравнения
///</summary>
complex<double> BoundaryValueProblem::NewtonMethod() {
	complex<double> a = IntegrandDenumerator(4);
	size_t iter = 0;
	while (abs(a) > eps && iter++ < 10) {
		alpha -= a;
		a = IntegrandDenumerator(4);
	}
	return alpha;
}



vector<complex<double>> BoundaryValueProblem::CauchyProblemSolutions(double y,
	size_t size) const {
	const complex<double> Mu = sqrt(Parameters::muConst);
	//const complex<double> kappa1 = kappa / Mu;
	const complex<double> mu = sqrt(alpha * alpha - kappa * kappa / Mu);
	const complex<double> cm = cosh(mu * y);
	const complex<double> sm = sinh(mu * y);
	vector<complex<double>> cauchyProblemSolution(size);
	cauchyProblemSolution[0] = sm / mu;
	cauchyProblemSolution[1] = cm * Parameters::muConst;
	const auto mu2 = mu * mu;
	if (size > 2) {
		cauchyProblemSolution[2] = alpha / mu2 * (cm * y - sm / mu);
		cauchyProblemSolution[3] = sm * alpha * y / mu * Parameters::muConst;
	}
	if (size > 4) {
		cauchyProblemSolution[4] = sm / mu / mu / mu * y * y * alpha * alpha -
			3.0 * cm / mu2 / mu2 * y * alpha * alpha + cm / mu2 * y +
			3.0 * sm / mu2 / mu2 / mu * alpha * alpha - sm / mu2 / mu;
		cauchyProblemSolution[5] = (cm * alpha * alpha * y * y / mu2 -
			sm * alpha * alpha * y / mu2 / mu + sm * y / mu) * Parameters::muConst;
	}
	return cauchyProblemSolution;
}

BoundaryValueProblem::BoundaryValueProblem()
{
}

BoundaryValueProblem::BoundaryValueProblem(double kappa, complex<double> alpha,
	double eps) :kappa(kappa), alpha(alpha),
	eps(eps)
{
}


BoundaryValueProblem::~BoundaryValueProblem() {}
