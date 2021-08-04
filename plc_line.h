#pragma once
#include <iostream>
#include <vector>

constexpr int DIMS = 3;

using dvec_t = std::vector<double>;
using bvec_t = std::vector<bool>;
using ivec_t = std::vector<int>;



void plc_line(const dvec_t& p1_in, const bvec_t& is_perf_in, dvec_t& p1_out, bvec_t& is_perf_out,
	double tol1, double tol2);

