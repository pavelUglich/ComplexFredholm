#include "Stabilizer.h"
#include<cmath>



void stabilizer::SquareRoot()
{
	Diagonal[0] = sqrt(Diagonal[0]);
	UpDiagonal[0] /= Diagonal[0];
	for (size_t i = 1; i < size - 1; i++) {
		Diagonal[i] = sqrt(Diagonal[i] - UpDiagonal[i - 1] * UpDiagonal[i - 1]);
		UpDiagonal[i] /= Diagonal[i];
	}
	Diagonal[size - 1] = sqrt(Diagonal[size - 1] - UpDiagonal[size - 2] * UpDiagonal[size - 2]);
}


stabilizer::stabilizer() :size(0)
{
}

stabilizer::stabilizer(size_t n, double step, double p, BoundaryCondition Left, BoundaryCondition Right) : size(n) {
	Diagonal = vector<double>(size);
	UpDiagonal = vector<double>(size - 1);
	double hStab = p / step / step;
	for (size_t i = 0; i < size; i++)  Diagonal[i] = 1 + 2 * hStab;
	for (size_t i = 0; i < size - 1; i++)  UpDiagonal[i] = -hStab;
	if (Left == Dirichle) Diagonal[0] = 1 + 3 * hStab;
	else Diagonal[0] = 1 + hStab;
	if (Right == Dirichle) Diagonal[size - 1] = 1 + 3 * hStab;
	else Diagonal[size - 1] = 1 + hStab;
	SquareRoot();
}

//~Stabilizer();

vector<double> stabilizer::diagonal() const { return Diagonal; }

vector<double> stabilizer::up_diagonal() const { return UpDiagonal; }

/*
Stabilizer::~Stabilizer()
{
}*/
