#pragma once
#include<vector>
#include<complex>
#include "Stabilizer.h"
//#include<boost/numeric/ublas/matrix.hpp>
using namespace std;

double norm(const vector<complex<double>> & v);
class matrix_system
{
	//size_t size;
	size_t _rows;
	size_t _columns;
	vector<vector<complex<double>>> _matrix;
	vector<complex<double>> RightPart;
	double step;
	vector<complex<double>> p1, p2;
	stabilizer _stabilizer;

	void multiply_ASinv();
	void del_col(size_t k);
	void del_row(size_t k);
	void multiply_transpose_au();
	void QPR();
	void multiply_rx();



public:
	//MatrixSystem();
	matrix_system(
		const vector<vector<complex<double>>> & A,
		const vector<complex<double>> & b,
		double step,
		double p = 1,
		BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann) :
		_matrix(A), RightPart(b), step(step) {
		_rows = _matrix.size();
		_columns = _matrix.front().size();

//		size = b.size();
		_stabilizer = stabilizer(_columns, step, p, left, right);
		multiply_ASinv();
		multiply_transpose_au();
		QPR();
		multiply_rx();
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

	vector<complex<double>> multiply_qtu(const vector<complex<double>> & v);
	void multiply_rtx(vector<complex<double>> &u);
	void multiply_sinv(vector<complex<double>> &u);
};

