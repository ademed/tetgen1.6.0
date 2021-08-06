#include "plc_line.h"
#include "templates_def.h"

double norm(const dvec_t& v1, const dvec_t& v2);

void plc_line(const dvec_t& p1_in, const bvec_t& is_perf_in, dvec_t& p1_out, bvec_t& is_perf_out,
	const double tol1, const double tol2)
{
	double C, S;												//C is straight line; S is curved edge
	dvec_t v1, v2;
	bool cond2;
	bool append_flag;
	size_t n_points = p1_in.size() / DIMS;
	if (n_points > 2)
	{
		//std::cout << p1_in.size() << '\n';
		for (size_t j = 0; j < DIMS; j++)
			p1_out.push_back(p1_in[j]);
		size_t i = 1;
		while (i < n_points)
		{
			if (is_perf_in[i] != is_perf_in[i - 1])
			{
				for (size_t j = i * DIMS; j < (i * DIMS + DIMS); j++)
					p1_out.push_back(p1_in[j]);
				is_perf_out.push_back(is_perf_in[i - 1]);
				i++;
			}
			v1 = slice(p1_in, i * DIMS, i * DIMS + DIMS - 1);
			v2 = slice(p1_out, p1_out.size() - DIMS, p1_out.size() - 1);
			S = norm(v1, v2);
			C = S;
			cond2 = (0.5 * sqrt(S * S - C * C)) < tol2;
			C = 0.;
			while ((i < n_points - 1) && (cond2) && (C < tol1))
			{
				append_flag = 0;
				if (is_perf_in[i] != is_perf_in[i - 1])
				{
					for (size_t j = i * DIMS; j < (i * DIMS + DIMS); j++)
						p1_out.push_back(p1_in[j]);
					is_perf_out.push_back(is_perf_in[i - 1]);
					append_flag = 1;
					i++;
					break;
				}
				i++;
				v1 = slice(p1_in, i * DIMS, i * DIMS + DIMS - 1);
				v2 = slice(p1_in, (i - 1) * DIMS, (i - 1) * DIMS + DIMS - 1);
				S += norm(v1, v2);
				v1 = slice(p1_in, i * DIMS, i * DIMS + DIMS - 1);
				v2 = slice(p1_out, p1_out.size() - DIMS, p1_out.size() - 1);
				C = norm(v1, v2);
				cond2 = (0.5 * sqrt(S * S - C * C)) < tol2;
			}
			if (i != (n_points - 1))
			{
				if (!append_flag)                  //what does the else of this statmnt falling through mean?
				{
					for (size_t j = (i - 1) * DIMS; j < ((i - 1) * DIMS + DIMS); j++)
						p1_out.push_back(p1_in[j]);
					is_perf_out.push_back(is_perf_in[i - 1]);
				}
			}
			else
			{
				if (cond2)
				{
					for (size_t j = i * DIMS; j < (i * DIMS + DIMS); j++)
						p1_out.push_back(p1_in[j]);
					is_perf_out.push_back(is_perf_in[i - 1]);
				}
				else
				{
					for (size_t j = (i - 1) * DIMS; j < ((i - 1) * DIMS + DIMS); j++)
						p1_out.push_back(p1_in[j]);
					is_perf_out.push_back(is_perf_in[i - 1]);
					for (size_t j = p1_in.size() - DIMS; j < p1_in.size(); j++)
						p1_out.push_back(p1_in[j]);
					is_perf_out.push_back(is_perf_in[i - 1]);
				}
				break;
			}

		}

	}

	else
	{
		p1_out = p1_in;
	}
}




double norm(const dvec_t& v1, const dvec_t& v2)
{
	//dvec_t v3;
	double sum{ 0 };

	for (size_t i = 0; i < v1.size(); i++) {
		sum += (v2[i] - v1[i]) * (v2[i] - v1[i]);               //change to pow
	}

	return sqrt(sum);
}