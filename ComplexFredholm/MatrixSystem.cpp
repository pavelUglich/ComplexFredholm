#include "MatrixSystem.h"
#include<vector>
#include<cassert>

double norm(const vector<complex<double>> & v) {
	double sum = 0;
	for (size_t i = 0; i < v.size(); i++)
		sum += v[i].real() * v[i].real() + v[i].imag()*v[i].imag();
	return sqrt(sum);
}

void normalize(vector<complex<double>> & v) {
	double NormV = norm(v);
	for (size_t i = 0; i < v.size(); i++)	v[i] /= NormV;
}

complex<double> innerprod(const vector<complex<double>> & a, 
	const vector<complex<double>> & b) {
	assert(a.size() == b.size());
	complex<double> sum = 0;
	for (size_t i = 0; i < a.size(); i++)
		sum += a[i] * conj(b[i]);
	return sum;
}

vector<complex<double>> matrix_system::multiply_qtu(
	const vector<complex<double>> & v)
{
	auto qtu = v;
	for (size_t i = 0; i < _columns; i++) {
		if (i > _rows)
		{
			break;
		}
		vector<complex<double>> a(_rows - i);
		for (size_t j = 0; j < _rows - i; j++) a[j] = _matrix[j + i][i];
		complex<double> sc = 0;
		for (size_t k = 0; k < _rows - i; k++)
			sc += conj(a[k]) * qtu[k + i];
		for (size_t j = i; j < qtu.size(); j++) qtu[j] -= 2.0 * a[j - i] * sc;
	}
	return qtu;
}

///<summary>
///Умножение справа на матрицу, обратную к матрице стабилизатора 
///</summary>
void matrix_system::multiply_ASinv()
{
	auto diagonal = _stabilizer.diagonal();
	auto up_diagonal = _stabilizer.up_diagonal();
	for (auto& i : _matrix) i[0] /= diagonal[0];
	for (size_t i = 1; i < _matrix[0].size(); i++)
	{
		for (auto& j : _matrix)
		{
			j[i] -= up_diagonal[i - 1] * j[i - 1];
			j[i] /= diagonal[i];
		}
	}
}

///<summary>
///Умножение слева на матрицу отражения
///</summary>
///<param name="k">k - номер столбца</param>
void matrix_system::del_col(size_t k)
{
	if (k >= _columns || k >= _rows)
	{
		return;
	}
	size_t l = _rows - k;
	vector<complex<double>> av(l);
	for (size_t i = 0; i < l; i++)
		av[i] = _matrix[i + k][k];//!!!
	av[0] -= norm(av);// *av[0] / abs(av[0]);
	normalize(av);
	//vv Поддиагональная часть столбца матрицы
	vector<complex<double>> vv(l);
	for (size_t i = 0; i < l; i++) vv[i] = _matrix[i + k][k];
	auto sc = innerprod(vv, av);
	auto pp = _matrix[k][k] - 2.0 * av[0] * sc;
	for (size_t i = k + 1; i < _columns; i++) {
		for (size_t j = 0; j < l; j++) vv[j] = _matrix[j + k][i];
		sc = innerprod(vv, av);
		for (size_t j = k; j < _rows; j++)
			_matrix[j][i] -= 2.0 * av[j - k] * sc;
	}
	for (size_t i = 0; i < l; i++) _matrix[i + k][k] = av[i];
	p1.push_back(pp);
}


///<summary>
///Умножение справа на матрицу отражения
///</summary>
///<param name="k">k - номер строки</param>
void matrix_system::del_row(size_t k)
{
	if (k >= _columns - 1 || k >= _rows)
	{
		return;
	}
	size_t l = _columns - k - 1;
	vector<complex<double>> av(l);
	for (size_t i = 0; i < l; i++) av[i] = _matrix[k][i + k + 1];
	av[0] -= norm(av);// *av[0] / abs(av[0]);
	normalize(av);
	vector<complex<double>> vv(l);
	for (size_t i = 0; i < l; i++) vv[i] = _matrix[k][i + k + 1];
	auto sc = innerprod(vv, av);
	auto pp = _matrix[k][k + 1] - 2.0 * av[0] * sc;
	for (size_t i = k + 1; i < _rows; i++)
	{
		for (size_t j = 0; j < l; j++) vv[j] = _matrix[i][j + k + 1];
		sc = innerprod(vv, av);
		for (size_t j = k + 1; j < _columns; j++)
			_matrix[i][j] -= 2.0 * av[j - k - 1] * sc;
	}
	for (size_t i = 0; i < l; i++) _matrix[k][i + k + 1] = av[i];
	p2.push_back(pp);
	//return pp;
}

///<summary>
///Проеобразование правой части, умножение на транспонированную матрицу СЛАУ;
///</summary>
void matrix_system::multiply_transpose_au()
{
	vector<complex<double>> v(_columns);
	for (size_t i = 0; i < _columns; i++) {
		v[i] = 0;
		for (size_t j = 0; j < _rows; j++)
			v[i] += conj(_matrix[j][i]) * RightPart[j];
	}
	RightPart = move(v);
}

///<summary>
///Сведение матрицы СЛАУ к двухдиагональному виду
///</summary>
void matrix_system::QPR()
{
	const auto size = std::min(_rows, _columns);
	p1.clear();
	p2.clear();
	for (size_t i = 0; i < size; i++)
	{
		del_col(i);
		del_row(i);
	}
}


///<summary>
///Умножение вектора правой части на матрицу, обратную к ортогональной матрице R;
///</summary>
void matrix_system::multiply_rx()
{
	for (size_t i = 0; i < _rows - 1; i++) {
		vector<complex<double>> av(_columns);
		for (size_t j = i + 1; j < _columns; j++)
			av[j] = _matrix[i][j];
		complex<double> sc = 0;
		for (size_t j = i + 1; j < _columns; j++)
			sc += av[j] * RightPart[j];
		for (size_t j = i + 1; j < _columns; j++)
			RightPart[j] -= 2.0 * conj(av[j]) * sc;
	}
}

///<summary>
///Умножение вектора u на матрицу, обратную к ортогональной матрице R;
///</summary>
///<param name="u">вектор правой части</param>
void matrix_system::multiply_rtx(vector<complex<double>> &u) {
	auto v = u;
	for (size_t i = 0; i < _rows; i++) {
		size_t l = _rows - i;
		if (_columns < l)
		{
			continue;
		}
		vector<complex<double>> a(_columns - l);
		for (size_t j = 0; j < i; j++)
			a[j] = _matrix[l - 1][j + l];
		complex<double> sc = 0;
		for (size_t j = 0; j < a.size(); j++)	
			sc += a[j] * v[j + l];
		for (size_t j = l; j < _columns; j++) 
			v[j] -= 2.0* conj(a[j - l]) * sc;
	}
	u = v;
}

///<summary>
///Умножение вектора u на матрицу, обратную к матрице стабилизатора;
///</summary>
///<param name="u">вектор правой части</param>
void matrix_system::multiply_sinv(vector<complex<double>>& u) {
	auto diagonal = _stabilizer.diagonal();
	auto up_diagonal = _stabilizer.up_diagonal();
	auto x = u;
	x[_columns - 1] = u[_columns - 1] / diagonal[_columns - 1];
	for (size_t i = 1; i < _columns; i++) {
		const size_t j = _columns - i - 1;
		x[j] = (u[j] - up_diagonal[j] * x[j + 1]) / diagonal[j];
	}
	u = x;
}
