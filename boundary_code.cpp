#include <iostream>
#include <iomanip>   //std::setprecision()
#include <string>
#include <vector>
#include <cstddef>      // std::size_t
#include <cmath>        
#include <algorithm>
#include <functional>
#include <fstream>  
#include <istream>  
#include <assert.h> 


template <typename T>
void printVec(const std::vector<T>& vec)
{
	for (size_t j = 0; j < vec.size(); j++)
		std::cout << vec[j] << ' ';
	std::cout << '\n';
}
//add two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	assert(a.size() == b.size());
	std::vector<T>result;
	result.reserve(a.size());
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::plus<T>());
	return result;
}

//subtract two vectors
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
	assert(a.size() == b.size());
	std::vector<T>result;
	result.reserve(a.size());
	std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<T>());
	return result;
}

//get length of line from two points
double length(std::vector<double> a)
{
	double sum{ 0 };
	//OpenMP parallel for??
	for (size_t i{ 0 }; i < a.size(); i++) {
		sum += (a[i] * a[i]);
	}
	return sqrt(sum);
}

std::vector<double> scalTimesVec(double scalar, const std::vector<double>& vec)
{
	//std::transform(vec.begin(), vec.end(), vec.begin(), [scalar](int& c) {return c * scalar; });
	std::vector<double> result;
	result.reserve(vec.size());
	for (size_t i{ 0 }; i < vec.size(); i++) {
		result.push_back(scalar * vec[i]);
	}
	return result;
}
void voidScalTimesVec(double scalar, std::vector<double>& vec)
{
	for (size_t i{ 0 }; i < vec.size(); i++) {
		vec[i] = scalar * vec[i];
	}
}

template <typename T>
void voidScalPlusVec(double scalar, std::vector<T>& vec)
{
	for (size_t i{ 0 }; i < vec.size(); i++) {
		vec[i] = scalar + vec[i];
	}
}

template <typename T>
std::vector<T> scalPlusVec(double scalar, const std::vector<T>& vec)
{
	std::vector<T> result(vec.size());
	//printVec(result);
	for (size_t i{ 0 }; i < vec.size(); i++) {
		result[i] = scalar + vec[i];               
	}
	return result;
}

template <typename T> 
void printPerLine(const std::vector<T> vec, int lineCount)
{
	for (size_t i{ 0 }; i < vec.size(); i++)
	{
		std::cout << std::right <<  vec[i] << "  ";
		if (!( (i+1) % lineCount) ) std::cout << '\n';
	}
}

template <>
void printPerLine(const std::vector<double> vec, int lineCount)
{
	for (size_t i{ 0 }; i < vec.size(); i++)
	{
		std::cout << std::fixed << std::setprecision(2) << vec[i] << "  ";
		if (!((i + 1) % lineCount)) std::cout << '\n';
	}
}
int main()
{
	//input parameters

	double tol = 0.4;
	std::vector<double>v0 = { 1, 0, 1 }; //intersect
	std::vector<double>v1 = { 0, 1, 0 };
	std::vector<double>v2 = { -1, 0, 0 };
	std::vector<double>v3 = { 0, 0, -1 };    //offset/displacement


	//constants
	constexpr int dim = 3;
	constexpr int quad = 4;

	unsigned long ndiv1 = ceil(length(v1) / tol) ;
	unsigned long ndiv2 = ceil(length(v2) / tol) ;
	unsigned long ndiv3 = ceil(length(v3) / tol) ;

	std::vector<double> v1_step, v2_step, v3_step;
	v1_step = scalTimesVec(1. / ndiv1, v1);
	v2_step = scalTimesVec(1. / ndiv2, v2);
	v3_step = scalTimesVec(1. / ndiv3, v3);

	//Note: Faces 1, 2 and 3 all share common vertex v0. Also, vectors v1, v2 v3 all share (and point outward) from common vertex v0
	std::vector<double> refined_face1, refined_face2, refined_face3, refined_face4, refined_face5, refined_face6;
	//v1, v2
	refined_face1.reserve((ndiv1+1)*(ndiv2+1));
	std::vector<double> start = v0;
		refined_face1.insert(std::end(refined_face1), std::begin(start), std::end(start));
	for (size_t i{ 1 }; i < (ndiv1 + 1)*(ndiv2 + 1); i++) 
	{
		std::vector<double> v = scalTimesVec(!(i % (ndiv1 + 1)) ? 1 : 0, v2_step);
		std::vector<double> tmp = start + scalTimesVec(i % (ndiv1 + 1), v1_step) + v ;
		start = start + v;
		refined_face1.insert(std::end(refined_face1), std::begin(tmp), std::end(tmp));
	}

	//v1, v3
	refined_face2.reserve((ndiv1 + 1) * (ndiv3 - 1));
	start = v0 + v3_step;
	refined_face2.insert(std::end(refined_face2), std::begin(start), std::end(start));
	for (size_t i{ 1 }; i < (ndiv1 + 1) * (ndiv3 - 1); i++)
	{
		std::vector<double> v = scalTimesVec(!(i % (ndiv1 + 1)) ? 1 : 0, v3_step);
		std::vector<double> tmp = start + scalTimesVec(i % (ndiv1 + 1), v1_step) + v;
		start = start + v;
		refined_face2.insert(std::end(refined_face2), std::begin(tmp), std::end(tmp));
	}

	//v2, v3
	refined_face3.reserve((ndiv2 - 1) * (ndiv3 - 1));
	start = v0 + v2_step + v3_step;
	refined_face3.insert(std::end(refined_face3), std::begin(start), std::end(start));
	for (size_t i{ 1 }; i < (ndiv2 - 1) * (ndiv3 - 1); i++)
	{
		std::vector<double> v = scalTimesVec(!(i % (ndiv2 - 1)) ? 1 : 0, v3_step);
		std::vector<double> tmp = start + scalTimesVec(i % (ndiv2 - 1), v2_step) + v;
		start = start + v;
		refined_face3.insert(std::end(refined_face3), std::begin(tmp), std::end(tmp));
	}

	//face1's opposite
	refined_face4 = refined_face1 ;
	for (size_t i{ 0 }; i < refined_face4.size(); i++) {
		refined_face4[i] +=  v3[i % (v3.size())];
	}

	//face2's opposite
	refined_face5 = refined_face2 ;
	for (size_t i{ 0 }; i < refined_face5.size(); i++) {
		refined_face5[i] +=  v2[i % (v2.size())];
	}

	//face3's opposite
	refined_face6 = refined_face3 ;
	for (size_t i{ 0 }; i < refined_face6.size(); i++) {
		refined_face6[i] +=  v1[i % (v1.size())];
	}

	//All points list
	std::vector<double> listPoints;
	listPoints.reserve(refined_face1.size() + refined_face2.size() + refined_face3.size() + refined_face4.size()
		+ refined_face5.size() + refined_face6.size());
	listPoints.insert(listPoints.end(), refined_face1.begin(), refined_face1.end());
	listPoints.insert(listPoints.end(), refined_face2.begin(), refined_face2.end());
	listPoints.insert(listPoints.end(), refined_face3.begin(), refined_face3.end());
	listPoints.insert(listPoints.end(), refined_face4.begin(), refined_face4.end());
	listPoints.insert(listPoints.end(), refined_face5.begin(), refined_face5.end());
	listPoints.insert(listPoints.end(), refined_face6.begin(), refined_face6.end());

	//create sub-facets from each of the six faces in (anticlockwise) order acc to listed points indices
	std::vector<int> facets1;
	facets1.reserve(8 * (ndiv1 * ndiv2 + ndiv1 * ndiv3 + ndiv2 * ndiv3));
	//refined_face1
	for (size_t i{ 0 }; i < (ndiv1 + 1) * ndiv2; i++)
	{
		if (!((i + 1) % (ndiv1 + 1))) {
			continue;
		}
		facets1.push_back(i);
		facets1.push_back(i + 1);
		facets1.push_back(i + ndiv1 + 2);
		facets1.push_back(i + ndiv1 + 1);
	}

	//facets4
	int offset4 = (ndiv1 + 1) * (ndiv2 + 1) + (ndiv1 + 1) * (ndiv3 - 1) + (ndiv2 - 1) * (ndiv3 - 1);
	std::vector<int> facets4 = scalPlusVec(offset4, facets1);

	//refined_face2
	std::vector<int> facets2;
	for (size_t i{ 0 }; i < (ndiv1 + 1) * (ndiv3); i++)
	{
		if (!((i + 1) % (ndiv1 + 1))) {
			continue;
		}
		facets2.push_back(i);
		facets2.push_back(i + 1);
		facets2.push_back(i + ndiv1 + 2);
		facets2.push_back(i + ndiv1 + 1);
	}

	double num = (ndiv1 + 1.) * (ndiv2 + 1.) - (ndiv1 + 1);
	voidScalPlusVec(num, facets2);
	for (size_t i = 0; i < quad * ndiv1; i += quad)
	{
		facets2[i] = facets1[i];
		facets2[i + 1] = facets1[i + 1];
	}

	for (size_t i = quad * (ndiv1) * (ndiv3 - 1) + (quad - 2), j = 1; i < facets2.size(); i += quad, j += quad)
	{
		facets2[i] = facets4[j];
		facets2[i + 1] = facets4[j - 1];
	}

	//refined_face5
	std::vector<int> facets5 = scalPlusVec(offset4, facets2);
	for (size_t i = 0, j = quad * (ndiv1) * (ndiv2 - 1) + (quad - 1); i < quad * ndiv1; i += quad, j += quad)
	{
		facets5[i] = facets1[j];
		facets5[i + 1] = facets1[j - 1];
	}

	for (size_t i = quad * (ndiv1) * (ndiv3 - 1) + (quad - 1), j = quad * (ndiv1) * (ndiv2 - 1) + (quad - 1);
		i < facets5.size(); i += quad, j += quad)
	{
		facets5[i] = facets4[j];
		facets5[i - 1] = facets4[j - 1]; 
	}


	//refined_face3
	std::vector<int> facets3(quad * ndiv2 * ndiv3);
	for (size_t i = 0, j = 0; i < quad * ndiv2; i += quad, j += (quad * ndiv1))
	{
		facets3[i] = facets1[j];
		facets3[i + 1] = facets1[j + quad - 1]; 
	}
	for (size_t i = facets3.size() - (quad * ndiv2) + quad - 1, j = 0; i < facets3.size(); i += quad, j += (quad * ndiv1))
	{
		facets3[i] = facets4[j];
		facets3[i - 1] = facets4[j + quad - 1];
	}
	for (size_t i = 0, j = 0; i < facets3.size(); i += (quad * ndiv2), j += (quad * ndiv1))
	{
		facets3[i] = facets2[j];
		facets3[i + quad - 1] = facets2[j + quad - 1]; 
	}

	for (size_t i = quad * ndiv2 - 3, j = 0; i < facets3.size(); i += (quad * ndiv2), j += (quad * ndiv1))
	{
		facets3[i] = facets5[j];
		facets3[i + 1] = facets5[j + quad - 1];
	}

	for (size_t i = 2, j = 0, k = (ndiv1 + 1) * (ndiv2 + 1) + (ndiv1 + 1) * (ndiv3 - 1);
		i < facets3.size() - quad * ndiv2; i += quad, k++)
	{
		if (!(i - (2 + quad * (ndiv2 - 1) + (quad * ndiv2 * j))))
		{
			j++; k--;
			continue;
		}
		facets3[i] = k;
		facets3[i + quad + 1] = k;
		facets3[i + ndiv2 * quad - 1] = k;
		facets3[i + ndiv2 * quad + 2] = k;
	}


	//refined_face6
	std::vector<int> facets6(facets3.size());
	for (size_t i = 0, j = quad * ndiv1 - 3; i < quad * ndiv2; i += quad, j += (quad * ndiv1))
	{
		facets6[i] = facets1[j];
		facets6[i + 1] = facets1[j + 1]; 
	}
	for (size_t i = facets6.size() - (quad * ndiv2) + quad - 1, j = quad * (ndiv1 - 1) + 1; i < facets6.size(); i += quad, j += (quad * ndiv1))
	{
		facets6[i] = facets4[j];
		facets6[i - 1] = facets4[j + 1]; 
	}
	for (size_t i = 0, j = quad * (ndiv1 - 1) + 1; i < facets6.size(); i += (quad * ndiv2), j += (quad * ndiv1))
	{
		facets6[i] = facets2[j];
		facets6[i + quad - 1] = facets2[j + 1]; 
	}

	for (size_t i = quad * ndiv2 - 3, j = quad * (ndiv1 - 1) + 1; i < facets6.size(); i += (quad * ndiv2), j += (quad * ndiv1))
	{
		facets6[i] = facets5[j];
		facets6[i + 1] = facets5[j + 1]; 
	}

	for (size_t i = 2, j = 0, k = 2 * (ndiv1 + 1) * (ndiv2 + 1) + 2 * (ndiv1 + 1) * (ndiv3 - 1) + (ndiv2 - 1) * (ndiv3 - 1);
		i < facets6.size() - quad * ndiv2; i += quad, k++)
	{
		if (!(i - (2 + quad * (ndiv2 - 1) + (quad * ndiv2 * j))))
		{
			j++; k--;
			continue;
		}
		facets6[i] = k;
		facets6[i + quad + 1] = k;
		facets6[i + ndiv2 * quad - 1] = k;
		facets6[i + ndiv2 * quad + 2] = k;
	}

	//All facets list
	std::vector<int> listFacets;
	listFacets.reserve(facets1.size() + facets2.size() + facets3.size() + facets4.size()
		+ facets5.size() + facets6.size());
	listFacets.insert(listFacets.end(), facets1.begin(), facets1.end());
	listFacets.insert(listFacets.end(), facets2.begin(), facets2.end());
	listFacets.insert(listFacets.end(), facets3.begin(), facets3.end());
	listFacets.insert(listFacets.end(), facets4.begin(), facets4.end());
	listFacets.insert(listFacets.end(), facets5.begin(), facets5.end());
	listFacets.insert(listFacets.end(), facets6.begin(), facets6.end());

	std::cout << "ndiv1 is " << ndiv1 << '\n';
	std::cout << "-------------------listPoints--------------------------------------" << '\n';
	//printPerLine(listPoints, 3);
	std::cout << '\n';
	std::cout << "-------------------listFacets--------------------------------------" << '\n';
	//printPerLine(listFacets, 4);
	int k = 0;
	std::ofstream ofile("bnd.smesh");
	ofile << listPoints.size()/dim << " " << dim << " " << 0 << " " << 0 << '\n';
	ofile << k << " ";
	for (size_t j = 0; j < listPoints.size(); j++) 
	{
		ofile << listPoints[j] << " ";
		if ((j + 1) % 3 == 0 && j < listPoints.size() - 1)
			ofile << '\n' << ++k << " ";
	}
	ofile << '\n' << listFacets.size()/quad << " " << 0 << '\n';
	ofile << quad << " ";
	for (size_t j = 0; j < listFacets.size(); j++) 
	{
		ofile << listFacets[j] << " ";
		if ((j + 1) % quad == 0  && j < listFacets.size() - 1)
			ofile << '\n' << quad << " ";
	}
	ofile << '\n' << 0 << '\n' << 0 ;
//std::cout << listPoints.size()/3 << '\n';

}
