#include<vector>
#include<complex>


using namespace std;

#pragma once


/**
 * \brief Краевая задача для однородного упругого слоя 
 * (решение выражается аналитически)
 */
class BoundaryValueProblem {
protected:

	/**
	 * \brief Решение вспомогательной задачи Коши
	 */
	vector<complex<double>> CauchyProblemSolution;

	/**
	 * \brief Погрешность расчетов
	 */
	double eps;

	/**
	 * \brief параметр преобразования Фурье
	 */
	complex<double> alpha;

	/**
	 * \brief Волновое число
	 */
	double kappa;
public:
	/**
	 * \brief Решение вспомогательной задачи Коши
	 * \param y 
	 * \param size 
	 * \return 
	 */
	virtual vector<complex<double>> CauchyProblemSolutions(double y, size_t size = 2) const;

	vector<vector<complex<double>>> InnerCauchyProblemSolution;	
	BoundaryValueProblem();

	/**
	 * \brief Конструктор 
	 * \param kappa волновое число 
	 * \param alpha параметр преобразования Фурье
	 * \param eps погрешность вычислений
	 */
	BoundaryValueProblem(double kappa, complex<double> alpha, double eps = 0.1e-5);
	/**
	 * \brief деструктор
	 */
	virtual ~BoundaryValueProblem();
	complex<double> NewtonMethod();
	vector<complex<double>> GetCauchyProblemSolution() const { return CauchyProblemSolution; };
	/**
	 * \brief Знаменатель подынтегрального выражения
	 * \param size размер системы дифференциальных уравнений
	 * \return знаменатель
	 */
	virtual complex<double> IntegrandDenumerator(size_t size);
	virtual complex<double> waveField() const;
	virtual complex<double> perturb(double y) const;
	void SetAlpha(complex<double> _alpha);
	virtual void CreateTheInnerSolution(size_t points, size_t size = 2);
	double getKappa() const;
};

