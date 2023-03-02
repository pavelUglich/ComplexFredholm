#pragma once
#include<vector>
#include<complex>
#include "Stabilizer.h"
//#include<boost/numeric/ublas/matrix.hpp>
using namespace std;

double norm(const vector<complex<double>> & v);
class matrix_system
{
	size_t size;
	vector<vector<complex<double>>> Matrix;
	vector<complex<double>> RightPart;
	double step;
	vector<complex<double>> p1, p2;
	Stabilizer stabilizer;

	void multiply_ASinv();
	complex<double> DelCol(size_t k);
	complex<double> DelRow(size_t k);
	void MultiplyTransposeAu();
	void QPR();
	void multiply_Rx();



public:
	//MatrixSystem();
	matrix_system(
		const vector<vector<complex<double>>> & A,
		const vector<complex<double>> & b,
		double step,
		double p = 1,
		BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann) :
		Matrix(A), RightPart(b), step(step) {
		size = b.size();
		stabilizer = Stabilizer(size, step, p, left, right);
		multiply_ASinv();
		MultiplyTransposeAu();
		QPR();
		multiply_Rx();
	};
	//~MatrixSystem();

	///<summary>
	///Диагональ двухдиагональной матрицы;
	///</summary>
	vector<complex<double>> Diagonal() const {
		return p1;
	};

	///<summary>
	///Наддиагональ двухдиагональной матрицы;
	///</summary>
	vector<complex<double>> UpDiagonal() const {
		return p2;
	};

	///<summary>
	///Правая часть СЛАУ;
	///</summary>
	vector<complex<double>> rightPart() const {
		return RightPart;
	};

	vector<complex<double>> MultiplyQtu(const vector<complex<double>> & v);
	void multiply_Rtx(vector<complex<double>> &u);
	void multiply_Sinv(vector<complex<double>> &u);
};

