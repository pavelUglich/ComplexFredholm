#pragma once
#include<vector>
#include "MatrixSystem.h"
#include "IterativeProcess.h"
//#include "Stabilizer.h"
using namespace std;

class voyevodin_method
{
	double AlphaInitialValue;//Начальное значение параметра регуляризации
	double step; //длина отрезка разбиения
	double eps;//Погрешность определения параметра регуляризации
	double h; //Погрешность оператора
	double delta;//Погрешность правой части
	vector<complex<double>> RightPart, p1, p2, Qtu;
	//size_t size;
	size_t _rows;
	size_t _columns;
	vector<complex<double>> Solution;

public:
	//voyevodin_method();

	/// <summary>
	/// Конструктор метода Воеводина
	/// </summary>
	/// <param name="matrix">Матрица СЛАУ;</param>
	/// <param name="rightpart">вектор правой части;</param>
	/// <param name="Step">Длина подотрезка;</param>
	/// <param name="Left">Краевое условие;</param>
	/// <param name="Right">Краевое условие;</param>
	/// <param name="p">Параметр стабилизатора;</param>
	/// <param name="alphaInitialValue">Начальное значение параметра регуляризации;</param>
	/// <param name="H">Погрешность оператора;</param>
	/// <param name="Delta">Погрешность правой части;</param>
	/// <param name="eps">Погрешность определения параметра регуляризации;</param>
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
		//1. Создаём систему и приводим её к двухдиагональному виду
		matrix_system ms = { matrix, RightPart, step, p, Left, Right };
		//2. Запускаем итерационный процесс
		IterativeProcess iterativeProcess = {
			ms.Diagonal(),
			ms.UpDiagonal(),
			ms.rightPart(),
			ms.MultiplyQtu(RightPart),
			AlphaInitialValue, step, h,	delta,	eps };
		//3. Получаем решение
		Solution = iterativeProcess.solution();
		//4. Возвращаемя к изначальным неизвестнымю
		ms.multiply_Rtx(Solution);
		ms.multiply_Sinv(Solution);
	};


	~voyevodin_method();
	vector<complex<double>> solution() const {
		return Solution;
	};

};

