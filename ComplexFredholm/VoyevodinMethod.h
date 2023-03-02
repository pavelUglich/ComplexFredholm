#pragma once
#include<vector>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
//#include "Stabilizer.h"
using namespace std;

class voyevodin_method
{
	double AlphaInitialValue;//��������� �������� ��������� �������������
	double step; //����� ������� ���������
	double eps;//����������� ����������� ��������� �������������
	double h; //����������� ���������
	double delta;//����������� ������ �����
	vector<complex<double>> RightPart, p1, p2, Qtu;
	//size_t size;
	size_t _rows;
	size_t _columns;
	vector<complex<double>> Solution;

public:
	//voyevodin_method();

	/// <summary>
	/// ����������� ������ ���������
	/// </summary>
	/// <param name="matrix">������� ����;</param>
	/// <param name="rightpart">������ ������ �����;</param>
	/// <param name="Step">����� ����������;</param>
	/// <param name="Left">������� �������;</param>
	/// <param name="Right">������� �������;</param>
	/// <param name="p">�������� �������������;</param>
	/// <param name="alphaInitialValue">��������� �������� ��������� �������������;</param>
	/// <param name="H">����������� ���������;</param>
	/// <param name="Delta">����������� ������ �����;</param>
	/// <param name="eps">����������� ����������� ��������� �������������;</param>
	voyevodin_method(const vector<vector<complex<double>>> & matrix,
		const vector<complex<double>> & rightpart,
		double Step,
		BoundaryCondition Left = Neumann,
		BoundaryCondition Right = Neumann,
		double p = 1.0,
		double alphaInitialValue = 0.1e-5,
		double H = 0,
		double Delta = 0,
		double eps = 0.1e-11) :
		RightPart(rightpart), step(Step), AlphaInitialValue(alphaInitialValue), h(H), delta(Delta), eps(eps) {
		//1. ������ ������� � �������� � � ����������������� ����
		matrix_system ms = { matrix, RightPart, step, p, Left, Right };
		//2. ��������� ������������ �������
		IterativeProcess iterativeProcess = {
			ms.Diagonal(),
			ms.UpDiagonal(),
			ms.rightPart(),
			ms.MultiplyQtu(RightPart),
			AlphaInitialValue, step, h,	delta,	eps };
		//3. �������� �������
		Solution = iterativeProcess.solution();
		//4. ����������� � ����������� ������������
		ms.multiply_Rtx(Solution);
		ms.multiply_Sinv(Solution);
	};


	~voyevodin_method();
	vector<complex<double>> solution() const {
		return Solution;
	};

};

