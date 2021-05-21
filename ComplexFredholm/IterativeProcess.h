#pragma once
#include<vector>
#include<iostream>
#include<iterator>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
using namespace std;
class IterativeProcess
{
	vector<complex<double>> a;
	vector<complex<double>> b;
	vector<complex<double>> c;
	vector<complex<double>> p1;
	vector<complex<double>> p2;
	vector<complex<double>> Qtu;
	vector<complex<double>> RightPart;
	vector<complex<double>> Solution;
	double alpha;
	int Iterations;// = 50;
	size_t size;
	double step, h, delta, eps;

	///Перенести в итерационный процесс
	void tridiag();
	vector<complex<double>> Marching(const vector<complex<double>> & a);
	double Residual(vector<complex<double>> a, double alpha);
	void IterationsRun();

public:
	IterativeProcess(const vector<complex<double>> & p1, 
		const vector<complex<double>> & p2,
		const vector<complex<double>> & RightPart,
		const vector<complex<double>> & Qtu, double Alpha, double step, 
		double h, double delta,	double eps,	int iterations = 50) : p1(p1),	
		p2(p2),	RightPart(RightPart), Qtu(Qtu), alpha(Alpha), step(step), h(h),
		delta(delta), Iterations(iterations), eps(eps) {
		size = p1.size();
		tridiag();
		IterationsRun();
	};


	//~IterativeProcess();
	vector<complex<double>> solution() const {
		return Solution;
	}
};

