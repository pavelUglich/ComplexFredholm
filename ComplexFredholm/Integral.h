#pragma once
#include"BoundaryValueProblem.h"
#include <functional>
#include <utility>

const vector<double> xi{ -0.989400934991650 , -0.944575023073233, -0.865631202387832 ,
-0.755404408355003 , -0.617876244402644 , -0.458016777657227 , -0.281603550779259, -0.950125098376374e-1 ,
0.950125098376374e-1 ,0.281603550779259, 0.458016777657227, 0.617876244402644, 0.755404408355003,
0.865631202387832, 0.944575023073233, 0.989400934991650 };


const vector<double> eta{ 0.0271524594117540, 0.0622535239386475, 0.0951585116824939 ,
0.124628971255534 , 0.149595988816577 , 0.169156519395002 , 0.182603415044922, 0.189450610455068,
0.189450610455068 ,0.18260341504492, 0.16915651939500, 0.14959598881657, 0.124628971255534,
0.0951585116824939, 0.0622535239386475, 0.0271524594117540};

//complex<double> gauss(double a, double b, function<complex<double>(double)> contour, function<complex<double>(double)> derivative);

class Integral
{
	function<complex<double>()> integrand;
	BoundaryValueProblem *boundaryValueProblem;
	complex<double> gauss(double a, double b, double x,//,
		//function<complex<double>()> integrand = [this](){return boundaryValueProblem->waveField();},
		const function<complex<double>(double)>& contour = [](double y) {return complex<double>(y, 0); },
		const function<complex<double>(double)>& derivative = [](double y) {return complex<double>(1, 0); }) const;
public:
	complex<double> evaluate(double x) const;
	Integral();
	Integral(BoundaryValueProblem* boundaryValueProblem, function<complex<double>()> integrand);
	~Integral();
};

