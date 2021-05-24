#pragma once
#include<complex>
#include "BoundaryValueProblem.h"
using namespace std;


template<class BVP>
class Root {
	double kappa;
	complex<double> root;
	complex<double> Denumerator;
	complex<double> U;
	int N;
	vector<complex<double>> CauchyProblemSolution;

	vector<complex<double>> a01;
	vector<complex<double>> a11;

	vector<complex<double>> a02;
	vector<complex<double>> a12;

	vector<complex<double>> a03;
	vector<complex<double>> a13;

	complex<double> b1;
	complex<double> b2; 

public:


	vector<vector<complex<double>>> InnerCauchyProblemSolution;
	Root();
	Root(double kappa, complex<double> Root, size_t n);
	~Root();
	complex<double> Value() const { return root; }
	complex<double> u() const { return U; };
	complex<double> denumerator() const { return Denumerator; }
	void AddTheInnerSolution(size_t points);
	vector<complex<double>> GetCauchyProblemSolution() const { return CauchyProblemSolution; };


	vector<complex<double>> GetA01() const { return a01; }
	vector<complex<double>> GetA11() const { return a11; }

	vector<complex<double>> GetA02() const { return a02; }
	vector<complex<double>> GetA12() const { return a12; }

	vector<complex<double>> GetA03() const { return a03; }
	vector<complex<double>> GetA13() const { return a13; }

	complex<double> GetB1() const { return b1; }
	complex<double> GetB2() const { return b2; }

};


template<class BVP>
Root<BVP>::Root()
{
}

template<class BVP>
Root<BVP>::Root(double kappa, complex<double> Root_, size_t n) :kappa(kappa), N(n) {
	BVP boundaryBalueProblem = { kappa, Root_ };
 	root = boundaryBalueProblem.NewtonMethod();
	CauchyProblemSolution = boundaryBalueProblem.CauchyProblemSolutions(1, 6);
	Denumerator = CauchyProblemSolution[1];
	U = CauchyProblemSolution[0] / CauchyProblemSolution[3];
	b1 = CauchyProblemSolution[3];
	b2 = CauchyProblemSolution[5] / 2.0;
}


template<class BVP>
Root<BVP>::~Root() = default;

template<class BVP>
void Root<BVP>::AddTheInnerSolution(size_t points) {
	BoundaryValueProblem boundaryValueProblem = { kappa, root };
	boundaryValueProblem.CreateTheInnerSolution(points, 4);
	auto innerCauchyProblemSolution = boundaryValueProblem.InnerCauchyProblemSolution;
	for (size_t i = 0; i < points; i++)
	{
		a01.push_back(innerCauchyProblemSolution[i][1]);
		a11.push_back(innerCauchyProblemSolution[i][3]);
		a02.push_back(root * innerCauchyProblemSolution[i][0]);
		a12.push_back(innerCauchyProblemSolution[i][0] + root * innerCauchyProblemSolution[i][2]);
		a03.push_back(innerCauchyProblemSolution[i][0]);
		a13.push_back(innerCauchyProblemSolution[i][2]);
	}
}

