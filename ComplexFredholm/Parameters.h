#pragma once

#include<functional>
#include<complex>
#include<vector>
#include"Parameters.h"

using namespace std;

class Parameters {
public:
	static function<complex<double>(double)> Mu;
	static function<double(double)> Rho;
	static vector<complex<double>> mu;
	static vector<double> rho;
	static complex<double> muConst;
	static double rhoConst;

	/*
	Parameters() {
		Mu = [](double x) {return complex<double>(1.0, 0.0); };
		Rho = [](double x) {return 1.0; };
	}*/
};