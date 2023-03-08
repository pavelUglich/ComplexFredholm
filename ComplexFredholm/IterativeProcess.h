#pragma once
#include<vector>
#include<iostream>
#include<iterator>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
using namespace std;
class iterative_process
{
	vector<complex<double>> _a;
	vector<complex<double>> _b;
	vector<complex<double>> _c;
	vector<complex<double>> _p1;
	vector<complex<double>> _p2;
	vector<complex<double>> _qtu;
	vector<complex<double>> _right_part;
	vector<complex<double>> _solution;
	double _alpha;
	int _iterations;// = 50;
	size_t _size;
	double step, h, delta, eps;

	///Перенести в итерационный процесс
	void tridiag();
	vector<complex<double>> _marching(const vector<complex<double>> & a);
	double _residual(vector<complex<double>> a, double alpha);
	void iterations_run();

public:
	iterative_process(const vector<complex<double>> & p1, 
		const vector<complex<double>> & p2,
		const vector<complex<double>> & RightPart,
		const vector<complex<double>> & Qtu, double Alpha, double step, 
		double h, double delta,	double eps,	int iterations = 50);


	//~IterativeProcess();
	vector<complex<double>> solution() const {
		return _solution;
	}
};

