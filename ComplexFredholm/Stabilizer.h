#pragma once
#include<vector>
using namespace std;

enum BoundaryCondition { Dirichle, Neumann };


class stabilizer
{
	///Размер матрицы стабилизатора
	size_t size;
	///Диагональ
	vector<double> _diagonal;
	///Наддиагональ
	vector<double> _up_diagonal;
	/// Обработка симметричной трёхдиагональной матрицы по формулам метода квадратного корня
	void SquareRoot();

public:
	///конструкторы
	stabilizer();

	stabilizer(size_t n, double step, double p, BoundaryCondition Left = Neumann, BoundaryCondition Right = Neumann);;
	//~Stabilizer();
	vector<double> diagonal() const;
	vector<double> up_diagonal() const;
};


