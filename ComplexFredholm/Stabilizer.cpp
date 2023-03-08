#include "Stabilizer.h"
#include<cmath>



void stabilizer::SquareRoot()
{
	_diagonal[0] = sqrt(_diagonal[0]);
	_up_diagonal[0] /= _diagonal[0];
	for (size_t i = 1; i < size - 1; i++) {
		_diagonal[i] = sqrt(_diagonal[i] - _up_diagonal[i - 1] * _up_diagonal[i - 1]);
		_up_diagonal[i] /= _diagonal[i];
	}
	_diagonal[size - 1] = sqrt(_diagonal[size - 1] - _up_diagonal[size - 2] * _up_diagonal[size - 2]);
}


stabilizer::stabilizer() :size(0)
{
}

stabilizer::stabilizer(size_t n, double step, double p, BoundaryCondition Left, BoundaryCondition Right) : size(n) {
	_diagonal = vector<double>(size);
	_up_diagonal = vector<double>(size - 1);
	double hStab = p / step / step;
	for (size_t i = 0; i < size; i++)  _diagonal[i] = 1 + 2 * hStab;
	for (size_t i = 0; i < size - 1; i++)  _up_diagonal[i] = -hStab;
	if (Left == Dirichle) _diagonal[0] = 1 + 3 * hStab;
	else _diagonal[0] = 1 + hStab;
	if (Right == Dirichle) _diagonal[size - 1] = 1 + 3 * hStab;
	else _diagonal[size - 1] = 1 + hStab;
	SquareRoot();
}

//~Stabilizer();

vector<double> stabilizer::diagonal() const { return _diagonal; }

vector<double> stabilizer::up_diagonal() const { return _up_diagonal; }

/*
Stabilizer::~Stabilizer()
{
}*/
