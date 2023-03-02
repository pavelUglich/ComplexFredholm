#include "VoyevodinMethod.h"
#include "BoundaryValueProblem.h"
#include "BoundaryValueProblemSmooth.h"
#include "BoundaryValueProblemBrokenLine.h"
#include "SystemEvaluation.h"
#include "Parameters.h"
#include<iostream>
#include<complex>
#include<vector>
#include<clocale>
#include"Stabilizer.h"
#include <map>
#include "Integral.h"

using namespace std;

#include <fstream>


function<complex<double>(double)> Parameters::Mu;
function<double(double)> Parameters::Rho;
vector<complex<double>> Parameters::mu;
vector<double> Parameters::rho;
complex<double> Parameters::muConst;
double Parameters::rhoConst;

/**
 * \brief вывести на консоль значения вектора
 * \param points количество элементов
 * \param complexValuedVector вектор
 * \param fun лямбда-выражение
 * \param stream поток вывода
 */
void showVector(std::vector<std::complex<double>>& complexValuedVector,
	const function<double(complex<double>)>& fun, ostream& stream);

/**
 * \brief вывести на консоль значения функции
 * \param points количество точек
 * \param step шаг, расстояние между точками
 * \param complexValuedFunction фнкция
 * \param mapping отображение, применяемое к каждому значению функции
 * \param stream файл (или поток ввода-вывода)
 */
void showVector(size_t points, double step,
	const function<complex<double>(double)>& complexValuedFunction,
	const function<double(complex<double>)>& mapping, ostream& stream);


/**
 * \brief вывести в поток значения вектора
 * \param step шаг
 * \param complexValuedVector вектор
 * \param mapping отображение, которое применяется к каждому элементу вектора
 * \param stream поток вывода
 */
void showVector(double step,
	const std::vector<complex<double>>& complexValuedVector,
	const function<double(complex<double>)>& mapping, ostream& stream);


/**
 * \brief вывести в поток точное и восстановленное значения модуля сдвига
 * \param points количество точек
 * \param step шаг
 * \param fileName имя файла
 */
void plotTheSolution(size_t points, double step, const std::string& fileName);

/**
 * \brief Решить уравнение Фредгольма первого рода
 * \param Matrix аппроксимирующая матрица
 * \param rightPart вектор правой части
 * \param step шаг разбиентя
 * \return вектор решения
 */
vector<complex<double>> theFirstStep(const vector<vector<complex<double>>>& Matrix,
	const vector<complex<double>>& rightPart, double step);

/**
 * \brief пересчёт правой части уравнения
 * \param left вектор правой части, соответствующий текущей итерации
 * \param right вектор правой части, соответствующий наблюдаемому полю
 * \return новая правая часть
 */
vector<complex<double>> evaluateTheRightPart(
	const vector<complex<double>>& left, const vector<complex<double>>& right);

/**
 * \brief пересчёт модуля сдвига
 * \param solution вектор, добавляемый к текущему значению
 */
void reEvaluateParameters(const vector<complex<double>>& solution);

/**
 * \brief построить график волнового поля
 * \param step шаг по продольной координате
 * \param fieldHomogeneous
 * \param fieldCalculated
 * \param fieldObserved
 * \param fileName имя файла
 */
void plotTheWaveField(double step,
	const vector<complex<double>>& fieldHomogeneous,
	const vector<complex<double>>& fieldCalculated,
	const vector<complex<double>>& fieldObserved,
	const std::string& fileName);


void plotTheWaveField(double step,
	const std::map<string, vector<complex<double>>>& vectors,
	const std::string& fileName);

/**
 * \brief заполнить вектор значениями функции
 * \param points количество точек
 * \param step шаг по продольной координате
 * \param function функция
 * \return вектор
 */
vector<complex<double>> fillUpTheVector(size_t points, double step, 
	const function<complex<double>(double)>& function);

void setParameters(const size_t points, const double kappa, const double tau,
	const double g, const double h, const function<complex<double>(double)>& gx,
	const function<complex<double>(double)>& hx, const function<double(double)>& rho)
{
	const complex<double> I(0, 1);
	Parameters::muConst = (h + I * tau * kappa * g) / (tau * kappa * I + 1.0);
	Parameters::rhoConst = 1;
	for (size_t i = 0; i < points; i++)
	{
		Parameters::mu.push_back(Parameters::muConst);
		Parameters::rho.push_back(1);
	}
	Parameters::Mu = [=](double t)
	{
		return (hx(t) + I * tau * kappa * gx(t)) / (tau * kappa * I + 1.0);
	};
	Parameters::Rho = rho;
}



map<double, vector<complex<double>>> dispersionalSet(
	size_t points, double theHighestFrequency, double h, double tau, double g, const function<complex<double>(double)>& gx,
	const function<complex<double>(double)>& hx);

void plotTheDispersionalCurves(
	map<double, vector<complex<double>>>& dispersionSet,
	const string& fileName);

vector<complex<double>> evaluateTheRightPartRho(const vector<complex<double>>& left, const vector<complex<double>>& right);

void reEvaluateParametersRho(const vector<complex<double>>& v);

void plotTheSolutionRho(size_t points, double step, const string& fileName);

int main()
{
	setlocale(LC_ALL, "Russian");

	const size_t size_kappa = 20;
	const size_t size_gamma = 20;
	const size_t size_s = 20; // количество столбцов
	const double h_kappa = 1.0 / size_kappa;
	const double h_gamma = 1.0 / size_gamma;
	const double h_s = 1.0 / size_s;

	std::vector<std::vector<double>> matrix(size_kappa);
	const std::function<double(double, double)> kernel = [](auto x, auto s) {return 1.0 / (1 + 10 * (x - s) * (x - s)); };
	const std::function<double(double)> exact_solution = [](auto s) {return sin(pi * s); };
	for (size_t i = 0; i < size_kappa; i++)
	{
		std::vector<double> row(size_s);
		const double x = (i + 0.5) * h_kappa;
		for (size_t ii = 0; ii < size_s; ii++)
		{
			const double s = (ii + 0.5) * h_s;
			row[ii] = h_s * kernel(x, s);
		}
		matrix[i] = row;
	}

	std::vector<double> exact(size_s);
	for (size_t i = 0; i < size_s; i++)
	{
		const double s = (i + 0.5) * h_s;
		exact[i] = exact_solution(s);
	}

	const auto right_part = matrix * exact;
	VoyevodinMethod voyevodin_method = { matrix, right_part, h_kappa, Dirichle, Dirichle };
	const auto solution = voyevodin_method.solution();


	/*
	const size_t points = 40;
	const double kappa = 1.0;
	const double tau = 0.1;
	const double g = 0.5;
	const double h = 1.0;
	const double rightBoundary = 0.5;

	///const auto step = rightBoundary / points;
	// дисперсионные кривые
/*	setParameters(points, kappa, tau, g, h, [](double t) {return 0.5 * exp(-t); }, [](double t) {return 1.0 - 0.2 * exp(t); }, [](double x) {return 1.0; });
	/*BoundaryValueProblem* boundary_value_problem = new BoundaryValueProblemSmooth(kappa, { 1, 0 });
	const auto integrand = [=]() { return boundary_value_problem->waveField(); };
	Integral* integral = new Integral(boundary_value_problem, integrand);
	
	vector<complex<double>> numerical;
	for (size_t i = 0; i < points; i++) {
		cout << (i + 0.5) * step << endl;
		numerical.push_back(integral->evaluate((i + 0.5) * step));
	}
	delete boundary_value_problem;*/
/*	const auto systemEvaluationSmooth = new SystemEvaluation<BoundaryValueProblemSmooth>(kappa, points);	
	const vector<complex<double>> fieldObserved = systemEvaluationSmooth->EvaluateTheRightPart(0, rightBoundary, points);
	//plotTheWaveField(step, { {"black", numerical}, {"blue", fieldObserved} }, "wavefield.txt");

	
	auto dispSet = dispersionalSet(points, 10, h, tau, g, [](double t) { return 0.5 * exp(-t); },
	                               [](double t) { return 1.0 + 0.2 * t; });
	plotTheDispersionalCurves(dispSet, "dispSet.txt");//*/
	//setParameters(points, kappa, tau, g, h, [](double t) {return 0.5*exp(-t); }, [](double t) {return 1.0 + 0.4 * sin(Pi*t); });
	// решение обратной задачи
	/*
	setParameters(points, kappa, tau, g, h, [](double t) {return 0.5 * exp(-t); },
		[](double t) {return 1 + 0.5 * sin(Pi * t); }, [](auto x) {return exp(-0.5*x); });

	const double step = 1.0 / points;

	// точное решение задачи
	auto exactSolution = fillUpTheVector(points, step, Parameters::Mu);

	// объект для построения системы и правой части для постоянного нулевого начального приближения
	SystemEvaluation<BoundaryValueProblem> systemEvaluation = { kappa, points };
	auto Matrix = systemEvaluation.EvaluateTheKernel(0, 1, 0, 1, points);
	const auto fieldHomogeneous = systemEvaluation.EvaluateTheRightPart(0, 1, points);

	// создаём новый объект для отыскания наблюдаемого поля перемещений и построения правой части
	const SystemEvaluation<BoundaryValueProblemSmooth> systemEvaluationSmooth = { kappa, points };
	const vector<complex<double>> fieldObserved = systemEvaluationSmooth.EvaluateTheRightPart(0, 1, points);
	vector<complex<double>> rightPart = evaluateTheRightPart(fieldHomogeneous, fieldObserved);
	cout << endl << Norma(rightPart) << endl;

	auto solution = theFirstStep(Matrix, rightPart, step);
	reEvaluateParameters(solution);
	//reEvaluateParametersRho(solution);

	SystemEvaluation<BoundaryValueProblemBrokenLine> systemEvaluationBroken = { systemEvaluationSmooth.getRoots(), kappa, points };
	Matrix = systemEvaluationBroken.EvaluateTheKernelRho(0, 1, 0, 1, points); 	
	auto fieldCalculated = systemEvaluationBroken.EvaluateTheRightPart(0, 1, points);
	rightPart = evaluateTheRightPartRho(fieldCalculated, fieldObserved);
	cout << endl << Norma(rightPart) << endl;
	solution = theFirstStep(Matrix, rightPart, step);
	reEvaluateParametersRho(solution);
	double nz;
	int iterations = 0;
	do
	{		
		systemEvaluationBroken = { systemEvaluationBroken.getRoots(), kappa, points };
		Matrix = systemEvaluationBroken.EvaluateTheKernel(0, 1, 0, 1, points);
		fieldCalculated = systemEvaluationBroken.EvaluateTheRightPart(0, 1, points);
		rightPart = evaluateTheRightPart(fieldCalculated, fieldObserved);
		cout << endl << Norma(rightPart) << endl;
		solution = theFirstStep(Matrix, rightPart, step);
		reEvaluateParameters(solution);
				
		systemEvaluationBroken = { systemEvaluationBroken.getRoots(), kappa, points };
		Matrix = systemEvaluationBroken.EvaluateTheKernelRho(0, 1, 0, 1, points);
		fieldCalculated = systemEvaluationBroken.EvaluateTheRightPart(0, 1, points);
		rightPart = evaluateTheRightPartRho(fieldCalculated, fieldObserved);
		solution = theFirstStep(Matrix, rightPart, step);
		reEvaluateParametersRho(solution);		
		nz = Norma(rightPart);
		cout << endl << nz << " " << ++iterations << endl;
	}
	while (nz > 0.1e-3 && iterations < 10);
	plotTheWaveField(step, fieldHomogeneous, fieldCalculated, fieldObserved, "wavefield.txt");
	plotTheSolution(points, step, "solution.txt");//
	plotTheSolutionRho(points, step, "solution1.txt");//*/
	system("pause");
	return 0;
}

void showVector(std::vector<std::complex<double>>& complexValuedVector,
	const function<double(complex<double>)>& fun, ostream& stream)
{
	const auto points = complexValuedVector.size();
	for (size_t i = 0; i < points; i++)
	{
		stream << "[" << (i + 0.5) / points << "," << fun(complexValuedVector[i]) << "],";
	}
	cout << endl;
}

void showVector(size_t points, double step, 
	const function<complex<double>(double)>& complexValuedFunction, 
	const function<double(complex<double>)>& mapping, ostream& stream)
{
	for (size_t i = 0; i < points; i++)
	{
		const auto x = (i + 0.5) * step;
		stream << "(" << x << ", " << mapping(complexValuedFunction(x)) << ") ";
	}
}

void showVector(double step,
	const std::vector<complex<double>>& complexValuedVector,
	const function<double(complex<double>)>& mapping, ostream& stream)
{
	for (size_t i = 0; i < complexValuedVector.size(); i++)
	{
		const auto x = (i + 0.5) * step;
		stream << "(" << x << ", " << mapping(complexValuedVector[i]) << ") ";
	}
}

void showVector(double step,
                const std::vector<double>& vec,
                ostream& stream)
{
	for (size_t i = 0; i < vec.size(); i++)
	{
		const auto x = (i + 0.5) * step;
		stream << "(" << x << ", " << vec[i] << ") ";
	}
}


void plotTheSolution(size_t points, double step, const std::string& fileName)
{
	ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	stream << "		\\addplot[smooth, mark = *, blue] plot coordinates{\n";
	showVector(points, step, Parameters::Mu, [](auto x) { return x.real(); }, stream);
	stream << "					};\n";
	stream << "		\\addplot[line width = 0.25mm, smooth, black] plot coordinates{\n";
	showVector(step, Parameters::mu, [](auto x) { return x.real(); }, stream);
	stream << "					};\n";
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	stream << "		\\addplot[smooth, mark = *, red] plot coordinates{\n";
	showVector(step, Parameters::mu, [](auto x) { return x.imag(); }, stream);
	stream << "					};\n";
	stream << "		\\addplot[line width = 0.25mm, smooth, black] plot coordinates{\n";
	showVector(points, step, Parameters::Mu, [](auto x) { return x.imag(); }, stream);
	stream << "					};\n";
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}

vector<complex<double>> theFirstStep(
	const vector<vector<complex<double>>>& Matrix, 
	const vector<complex<double>>& rightPart, double step)
{
	const VoyevodinMethod voyevodinMethod(Matrix, rightPart, step, 
		Dirichle, Dirichle, 2, 0.1e-3, 0, 0, 
		0.1e-5);
	return voyevodinMethod.solution();
}

vector<complex<double>> evaluateTheRightPart(
	const vector<complex<double>>& left, 
	const vector<complex<double>>& right)
{
	assert(left.size() == right.size());
	vector<complex<double>> result(left.size());
	for (size_t i = 0; i < left.size(); i++)
	{
		result[i] = right[i] - left[i];
	}
	return result;
}

void reEvaluateParameters(const vector<complex<double>>& solution)
{
	assert(Parameters::mu.size() == solution.size());
	for (size_t i = 0; i < solution.size(); i++)
	{
		Parameters::mu[i] += solution[i];
	}
}

void addTheCurve(double step, const vector<complex<double>>& fieldHomogeneous,
	ofstream& stream, const std::string& color)
{
	stream << "	\\addplot[line width = 0.25mm, smooth, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(step, fieldHomogeneous, [](auto x) { return x.real(); }, stream);
	stream << "					};\n";
	stream << "		\\addplot[smooth, mark = *, ";
	stream << color;
	stream << "] plot coordinates{\n";
	showVector(step, fieldHomogeneous, [](auto x) { return x.imag(); }, stream);
	stream << "					};\n";
}

void plotTheWaveField(double step,
	const vector<complex<double>>& fieldHomogeneous,
	const vector<complex<double>>& fieldCalculated,
	const vector<complex<double>>& fieldObserved, const std::string& fileName)
{
	ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	addTheCurve(step, fieldHomogeneous, stream, "black");
	addTheCurve(step, fieldCalculated, stream, "blue");
	addTheCurve(step, fieldObserved, stream, "red");
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}

void plotTheWaveField(double step, const std::map<string, 
	vector<complex<double>>>& vectors, const std::string& fileName)
{
	ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	for(const auto& item: vectors)
	{
		addTheCurve(step, item.second, stream, item.first);
	}
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();
}


vector<complex<double>> fillUpTheVector(const size_t points, double step,
	const function<complex<double>(double)>& function)
{
	vector<complex<double>> exactSolution;
	for (size_t i = 0; i < points; i++)
	{
		const auto x = (i + 0.5) * step;
		exactSolution.emplace_back(function(x));
	}
	return exactSolution;
}

map<double, vector<complex<double>>> dispersionalSet(size_t points, 
	double theHighestFrequency, double h, double tau, double g, 
	const function<complex<double>(double)>& gx,
	const function<complex<double>(double)>& hx)
{
	map<double, vector<complex<double>>> result;
	double kappa = 0.1;
	setParameters(points, kappa, tau, g, h, gx, hx, [] (double x) {return 1.0; });
	auto system_evaluation = new SystemEvaluation<BoundaryValueProblemSmooth>
		(kappa, points);
	auto roots = system_evaluation->getRoots();
	while (kappa < theHighestFrequency)
	{
		const double step = 0.1;
		cout << kappa << endl;
		result[kappa] = roots;
		kappa += step;
		setParameters(points, kappa, tau, g, h, gx, hx, []( auto x ) {return 1; });
		delete system_evaluation;
		system_evaluation = new SystemEvaluation<BoundaryValueProblemSmooth>(roots, kappa, points);
		roots = system_evaluation->getRoots();
	}
	delete system_evaluation;
	return result;
}

void plotTheDispersionalCurves(
	map<double, vector<complex<double>>>& dispersionSet, 
	const string& fileName)
{
	ofstream os(fileName);
	os << "\\begin{tikzpicture}[scale=1.5]\n";
	os << "\\begin{axis}[grid]\n";
	const auto numberOfCurves = dispersionSet.begin()->second.size();
	for (auto& item : dispersionSet) {
		sort(item.second.begin(), item.second.end(), 
			[](auto x, auto y) { return x.real() < y.real(); });
	}
	const auto maxFreq = dispersionSet.rbegin()->first;
	for (size_t i = 0; i < numberOfCurves; i++) {
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto x : dispersionSet) {
			os << "(" << x.second[i].real() << ", " << x.first << ") ";
		}
		os << "};\n";
		os << "\\addplot[smooth, black] plot coordinates{\n";
		for (auto x : dispersionSet) {
			const auto value = x.second[i].imag();
			if (value < maxFreq) {
				os << "(" << -value << ", " << x.first << ") ";
			}
		}
		os << "};\n";
	}
	os << "\\end{axis}\n";
	os << "\\end{tikzpicture}\n";
	os.close();
}

vector<complex<double>> evaluateTheRightPartRho(const vector<complex<double>>& left, const vector<complex<double>>& right)
{
	assert(left.size() == right.size());
	vector<complex<double>> result(left.size());
	for (size_t i = 0; i < left.size(); i++)
	{
		result[i] = { (right[i] - left[i]).real(),0 };
	}
	return result;
}

void reEvaluateParametersRho(const vector<complex<double>>& v)
{
	assert(Parameters::rho.size() == v.size());
	for (size_t i = 0; i < v.size(); i++)
	{
		Parameters::rho[i] += v[i].real();
	}
}

void plotTheSolutionRho(size_t points, double step, const string& fileName)
{
	ofstream stream(fileName);
	stream << "\\begin{tikzpicture}[scale=1.5]\n";
	stream << "\\begin{axis}[grid]\n";
	stream << "		\\addplot[smooth, mark = *, blue] plot coordinates{\n";
	showVector(points, step, Parameters::Rho, [](auto x) { return x.real(); }, stream);
	stream << "					};\n";
	stream << "		\\addplot[line width = 0.25mm, smooth, black] plot coordinates{\n";
	showVector(step, Parameters::rho, stream);
	stream << "					};\n";
	stream << "\\end{axis}\n";
	stream << "\\end{tikzpicture}\n";
	stream.close();

}
