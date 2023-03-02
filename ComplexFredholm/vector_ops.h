#pragma once
#include <vector>

template<class T>
std::vector<T> operator*(const std::vector<std::vector<T>>& matrix,
	const std::vector<double>& vector)
{
	std::vector<T> result(matrix.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		for (size_t ii = 0; ii < vector.size(); ii++)
		{
			result[i] += matrix[i][ii] * vector[ii];
		}
	}
	return result;
}
