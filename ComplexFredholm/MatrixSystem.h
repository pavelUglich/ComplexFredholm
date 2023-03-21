#pragma once
#include<vector>
#include<complex>
#include "Stabilizer.h"


double norm(const vector<complex<double>> & v);

class matrix_system
{
	size_t _rows;
	size_t _columns;
	vector<vector<complex<double>>> _matrix;
	vector<complex<double>> _right_part;
	double _step;
	vector<complex<double>> _p1;
	vector<complex<double>> _p2;
	stabilizer _stabilizer;

	void multiply_ASinv();
	void del_col(size_t k);
	void del_row(size_t k);
	void multiply_transpose_au();
	void qpr();
	void multiply_rx();



public:
	matrix_system(
		const vector<vector<complex<double>>> & A,
		const vector<complex<double>> & b,
		double step,
		double p = 1,
		BoundaryCondition left = Neumann,
		BoundaryCondition right = Neumann);

	///<summary>
	///Диагональ двухдиагональной матрицы;
	///</summary>
	vector<complex<double>> Diagonal() const {
		return _p1;
	};

	///<summary>
	///Наддиагональ двухдиагональной матрицы;
	///</summary>
	vector<complex<double>> UpDiagonal() const {
		return _p2;
	};

	///<summary>
	///Правая часть СЛАУ;
	///</summary>
	vector<complex<double>> rightPart() const {
		return _right_part;
	};

	vector<complex<double>> multiply_qtu(const vector<complex<double>> & v);
	void multiply_rtx(vector<complex<double>> &u);
	void multiply_sinv(vector<complex<double>> &u);
};

