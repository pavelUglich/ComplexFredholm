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
 * \brief ������� �� ������� �������� �������
 * \param points ���������� ���������
 * \param complexValuedVector ������
 * \param fun ������-���������
 * \param stream ����� ������
 */
void showVector(std::vector<std::complex<double>>& complexValuedVector,
	const function<double(complex<double>)>& fun, ostream& stream);

/**
 * \brief ������� �� ������� �������� �������
 * \param points ���������� �����
 * \param step ���, ���������� ����� �������
 * \param complexValuedFunction ������
 * \param mapping �����������, ����������� � ������� �������� �������
 * \param stream ���� (��� ����� �����-������)
 */
void showVector(size_t points, double step,
	const function<complex<double>(double)>& complexValuedFunction,
	const function<double(complex<double>)>& mapping, ostream& stream);


/**
 * \brief ������� � ����� �������� �������
 * \param step ���
 * \param complexValuedVector ������
 * \param mapping �����������, ������� ����������� � ������� �������� �������
 * \param stream ����� ������
 */
void showVector(double step,
	const std::vector<complex<double>>& complexValuedVector,
	const function<double(complex<double>)>& mapping, ostream& stream);


/**
 * \brief ������� � ����� ������ � ��������������� �������� ������ ������
 * \param points ���������� �����
 * \param step ���
 * \param fileName ��� �����
 */
void plotTheSolution(size_t points, double step, const std::string& fileName);

/**
 * \brief ������ ��������� ���������� ������� ����
 * \param Matrix ���������������� �������
 * \param rightPart ������ ������ �����
 * \param step ��� ���������
 * \return ������ �������
 */
vector<complex<double>> theFirstStep(const vector<vector<complex<double>>>& Matrix,
	const vector<complex<double>>& rightPart, double step);

/**
 * \brief �������� ������ ����� ���������
 * \param left ������ ������ �����, ��������������� ������� ��������
 * \param right ������ ������ �����, ��������������� ������������ ����
 * \return ����� ������ �����
 */
vector<complex<double>> evaluateTheRightPart(
	const vector<complex<double>>& left, const vector<complex<double>>& right);

/**
 * \brief �������� ������ ������
 * \param solution ������, ����������� � �������� ��������
 */
void reEvaluateParameters(const vector<complex<double>>& solution);

/**
 * \brief ��������� ������ ��������� ����
 * \param step ��� �� ���������� ����������
 * \param fieldHomogeneous
 * \param fieldCalculated
 * \param fieldObserved
 * \param fileName ��� �����
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
 * \brief ��������� ������ ���������� �������
 * \param points ���������� �����
 * \param step ��� �� ���������� ����������
 * \param function �������
 * \return ������
 */
vector<complex<double>> fillUpTheVector(const size_t points, double step,
	const function<complex<double>(double)>& function);

void setParameters(const size_t points, const double kappa, const double tau,
	const double g, const double h, const function<complex<double>(double)>& gx,
	const function<complex<double>(double)>& hx)
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
	Parameters::Rho = [](double t)
	{
		return 1;
	};
}



map<double, vector<complex<double>>> dispersionalSet(
	size_t points, double theHighestFrequency, double h, double tau, double g, const function<complex<double>(double)>& gx,
	const function<complex<double>(double)>& hx);

void plotTheDispersionalCurves(
	map<double, vector<complex<double>>>& dispersionSet,
	const string& fileName);

int main()
{
	setlocale(LC_ALL, "Russian");
	const size_t points = 40;
	const double kappa = 2.0;
	const double tau = 0.0;
	const double g = 0.5;
	const double h = 1.0;
	const double rightBoundary = 0.5;

	/*const auto step = rightBoundary / points;
	setParameters(points, kappa, tau, g, h, [](double t) {return 0.5 * exp(-t); }, [](double t) {return 1.0 - 0.2 * exp(t); });
	BoundaryValueProblem* boundary_value_problem = new BoundaryValueProblemSmooth(kappa, { 1, 0 });
	const auto integrand = [=]() { return boundary_value_problem->waveField(); };
	Integral* integral = new Integral(boundary_value_problem, integrand);
	
	vector<complex<double>> numerical;
	for (size_t i = 0; i < points; i++) {
		cout << (i + 0.5) * step << endl;
		numerical.push_back(integral->evaluate((i + 0.5) * step));
	}
	delete boundary_value_problem;
	const auto systemEvaluationSmooth = new SystemEvaluation<BoundaryValueProblemSmooth>(kappa, points);	
	const vector<complex<double>> fieldObserved = systemEvaluationSmooth->EvaluateTheRightPart(0, rightBoundary, points);
	plotTheWaveField(step, { {"black", numerical}, {"blue", fieldObserved} }, "wavefield.txt");

	
	auto dispSet = dispersionalSet(points, 10, h, tau, g, [](double t) { return 0.5 * exp(-t); },
	                               [](double t) { return 1.0 + 0.2 * t; });
	plotTheDispersionalCurves(dispSet, "dispSet.txt");*/
	//setParameters(points, kappa, tau, g, h, [](double t) {return 0.5*exp(-t); }, [](double t) {return 1.0 + 0.4 * sin(Pi*t); });
	setParameters(points, kappa, tau, g, h, [](double t) {return 0.5 * exp(-t); }, [](double t) {return 1 + 0.5 * sin(Pi * t); });

	const double step = 1.0 / points;
	auto exactSolution = fillUpTheVector(points, step, Parameters::Mu);

	auto systemEvaluation = new SystemEvaluation<BoundaryValueProblem>(kappa, points);
	auto Matrix = systemEvaluation->EvaluateTheKernel(0, 1, 0, 1, points);
	const auto fieldHomogeneous = systemEvaluation->EvaluateTheRightPart(0, 1, points);

	SystemEvaluation<BoundaryValueProblemSmooth>* systemEvaluationSmooth = new SystemEvaluation<
		BoundaryValueProblemSmooth>(systemEvaluation->getRoots(), kappa, points);
	const vector<complex<double>> fieldObserved = systemEvaluationSmooth->EvaluateTheRightPart(0, 1, points);
	delete systemEvaluation;

	vector<complex<double>> rightPart = evaluateTheRightPart(fieldHomogeneous, fieldObserved);
	auto solution = theFirstStep(Matrix, rightPart, step);
	reEvaluateParameters(solution);
	SystemEvaluation<BoundaryValueProblemBrokenLine>* systemEvaluationBroken = new SystemEvaluation<
		BoundaryValueProblemBrokenLine>(systemEvaluationSmooth->getRoots(), kappa, points);
	delete systemEvaluationSmooth;
	Matrix = systemEvaluationBroken->EvaluateTheKernel(0, 1, 0, 1, points);
	auto fieldCalculated = systemEvaluationBroken->EvaluateTheRightPart(0, 1, points);
	rightPart = evaluateTheRightPart(fieldCalculated, fieldObserved);
	cout << endl << Norma(rightPart) << endl;
	double nz;
	int iterations = 0;
	do
	{
		solution = theFirstStep(Matrix, rightPart, step);
		reEvaluateParameters(solution);
		auto roots = systemEvaluationBroken->getRoots();
		delete systemEvaluationBroken;
		systemEvaluationBroken = new SystemEvaluation<BoundaryValueProblemBrokenLine>(
			roots, kappa, points);
		fieldCalculated = systemEvaluationBroken->EvaluateTheRightPart(0, 1, points);
		Matrix = systemEvaluationBroken->EvaluateTheKernel(0, 1, 0, 1, points);
		rightPart = evaluateTheRightPart(fieldCalculated, fieldObserved);
		nz = Norma(rightPart);
		cout << endl << nz << " " << ++iterations << endl;
	}
	while (nz > 0.05e-3 && iterations < 50);
	delete systemEvaluationBroken;
	plotTheWaveField(step, fieldHomogeneous, fieldCalculated, fieldObserved, "wavefield.txt");
	plotTheSolution(points, step, "solution.txt");
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
		0.1e-6);
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
	const double step = 0.1;
	setParameters(points, kappa, tau, g, h, gx, hx);
	auto system_evaluation = new SystemEvaluation<BoundaryValueProblemSmooth>
		(kappa, points);
	auto roots = system_evaluation->getRoots();
	while (kappa < theHighestFrequency)
	{
		cout << kappa << endl;
		result[kappa] = roots;
		kappa += step;
		setParameters(points, kappa, tau, g, h, gx, hx);
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
