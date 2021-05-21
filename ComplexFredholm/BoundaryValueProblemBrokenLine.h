#pragma once
#include "BoundaryValueProblemSmooth.h"
#include "Parameters.h"
#include<vector>


class BoundaryValueProblemBrokenLine :
	public BoundaryValueProblemSmooth{
	vector<complex<double>> CauchyProblem(double, double, vector<complex<double>> &) const override;


public:
	BoundaryValueProblemBrokenLine();
	BoundaryValueProblemBrokenLine(double kappa, complex<double> alpha, double eps = 0.1e-4);
	void CreateTheInnerSolution(size_t points, size_t size = 2) override;

	~BoundaryValueProblemBrokenLine();

};

