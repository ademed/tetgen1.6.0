#include <iostream>
#include <vector>
#include <iomanip>



template<typename T>
constexpr T PI{3.141592653589793238};

constexpr void linspace(int npoints, double first, double last, std::vector<double>& res){
    double dist = (last - first)/(npoints - 1); //get common difference between terms of arithmetic progression
    for (std::size_t i = 0; i < npoints; ++i)
    {
        res.push_back(first + i*dist);
    }   
}

template<typename T, typename S>
void CylinderPLC(T radius, S height, std::size_t refine = 100){ //refine defines the number of points 
    auto pi = PI<double>;
    std::vector<double> x,y,z,theta;
    theta.reserve(refine); x.reserve(refine); y.reserve(refine);z.reserve(refine);
    linspace(refine, 0.0, 2*PI<double>, theta); //populate theta with values btw 0 and 2pi using linear progression
    for(std::size_t i = 0; i < refine; ++i){
        x.push_back(radius*cos(theta[i]));
        y.push_back(radius*sin(theta[i]));
    }


}

int main(){


    #if 0  //testing linspace
    std::vector<double> test; test.reserve(10);
    linspace(20,1,10,test);
    for(auto m:test){
        std::cout << m << std::endl;
    }
    #endif

}