#include<vector>
#include<complex>


using namespace std;

#pragma once


/**
 * \brief ������� ������ ��� ����������� �������� ���� 
 * (������� ���������� ������������)
 */
class BoundaryValueProblem {
protected:

	/**
	 * \brief ������� ��������������� ������ ����
	 */
	vector<complex<double>> CauchyProblemSolution;

	/**
	 * \brief ����������� ��������
	 */
	double eps;

	/**
	 * \brief �������� �������������� �����
	 */
	complex<double> alpha;

	/**
	 * \brief �������� �����
	 */
	double kappa;
public:
	/**
	 * \brief ������� ��������������� ������ ����
	 * \param y 
	 * \param size 
	 * \return 
	 */
	virtual vector<complex<double>> CauchyProblemSolutions(double y, size_t size = 2) const;

	vector<vector<complex<double>>> InnerCauchyProblemSolution;	
	BoundaryValueProblem();

	/**
	 * \brief ����������� 
	 * \param kappa �������� ����� 
	 * \param alpha �������� �������������� �����
	 * \param eps ����������� ����������
	 */
	BoundaryValueProblem(double kappa, complex<double> alpha, double eps = 0.1e-5);
	/**
	 * \brief ����������
	 */
	virtual ~BoundaryValueProblem();
	complex<double> NewtonMethod();
	vector<complex<double>> GetCauchyProblemSolution() const { return CauchyProblemSolution; };
	/**
	 * \brief ����������� ���������������� ���������
	 * \param size ������ ������� ���������������� ���������
	 * \return �����������
	 */
	virtual complex<double> IntegrandDenumerator(size_t size);
	virtual complex<double> waveField() const;
	virtual complex<double> perturb(double y) const;
	void SetAlpha(complex<double> _alpha);
	virtual void CreateTheInnerSolution(size_t points, size_t size = 2);
	double getKappa() const;
};

