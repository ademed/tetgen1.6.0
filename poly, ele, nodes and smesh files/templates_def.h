#pragma once
#include <iostream>
#include <vector>

template<typename T>
//This slice INCLUDES the start and END points
std::vector<T> slice(std::vector<T> const& v, int m, int n)
{
	auto first = v.cbegin() + m;
	auto last = v.cbegin() + n + 1;

	std::vector<T> vec(first, last);
	return vec;
}


template <typename T>
void print_vec(const std::vector<T>& vec)
{
	for (size_t j = 0; j < vec.size(); j++)
		std::cout << vec[j] << ' ';
	std::cout << '\n';
}

template <typename T>
void print_mat(const std::vector<T>& vec)
{
	for (size_t j = 0; j < vec.size(); j++)
	{

		std::cout << vec[j] << ' ';
		if ( (j +1)% 3 == 0)
			std::cout << '\n';

	}
}