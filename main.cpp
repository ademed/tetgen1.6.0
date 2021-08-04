#include <iostream>
#include <string>
#include <vector>
#include <cstddef>      // std::size_t
#include <cmath>        // std::sqrt(double)
#include <fstream>  
#include <istream>  
#include "plc_line.h"
#include "templates_def.h"


int main()
{
	std::ifstream wfile("input.csv");

	dvec_t p1_in;
	bvec_t is_perf_in;
	std::string x_coord, y_coord, z_coord, perf;
	int i = 0;
	while (getline(wfile, x_coord, ','))
	{
		p1_in.push_back(stof(x_coord));
		getline(wfile, y_coord, ',');
		p1_in.push_back(stof(y_coord));
		getline(wfile, z_coord, ',');
		p1_in.push_back(stof(z_coord));
		if (i == 0) {
			std::string line;
			getline(wfile, line);
		}
		else {
			getline(wfile, perf, '\n');
			is_perf_in.push_back(stof(perf));
		}
		i++;
	}


	print_vec(p1_in);
	print_vec(is_perf_in);

	bvec_t is_perf_out;
	dvec_t p1_out; //Try reserving before pushing back...don't copy from stack to heap? Chern
	p1_out.reserve(p1_in.size());
	is_perf_out.reserve(is_perf_in.size());
	double tol1 = 4000.0;
	double tol2 = 120.0;
	plc_line(p1_in, is_perf_in, p1_out, is_perf_out, tol1, tol2);

	print_vec(p1_out);
	std::cout << "perforation values" << '\n';
	print_vec(is_perf_out);


	std::ofstream ofile("output.csv");
	for (size_t j = 0; j < p1_out.size(); j++) {
		ofile << p1_out[j] << ",";
		if ((j + 1) % 3 == 0)
			ofile << '\n';
	}
}