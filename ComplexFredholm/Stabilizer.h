#pragma once
#include<vector>
using namespace std;

enum BoundaryCondition { Dirichle, Neumann };


class stabilizer
{
	///������ ������� �������������
	size_t size;
	///���������
	vector<double> _diagonal;
	///������������
	vector<double> _up_diagonal;
	/// ��������� ������������ ��������������� ������� �� �������� ������ ����������� �����
	void SquareRoot();

public:
	///������������
	stabilizer();

	stabilizer(size_t n, double step, double p, BoundaryCondition Left = Neumann, BoundaryCondition Right = Neumann);;
	//~Stabilizer();
	vector<double> diagonal() const;
	vector<double> up_diagonal() const;
};


