#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>



template<typename T>
constexpr T PI{3.141592653589793238};

constexpr void linspace(int npoints, double first, double last, std::vector<double>& res){
    double dist = (last - first)/(npoints - 1); //get common difference between terms of arithmetic progression
    for (std::size_t i = 0; i < npoints; ++i)
    {
        res.push_back(first + i*dist);
    }   
}

template<typename T>
void write_smesh(T const& x_top, T const& y_top, T const& x_bot, T const& y_bot,
                std::vector<int> const& index_top,std::vector<int> const& index_bot, double h, std::size_t refine){
    
    std::ofstream ofile("cylinder.smesh"); 
    
    //write nodes
    std::size_t N = 2*refine;
    ofile << N << " " << 3 <<" " << 0 << " " << 0 <<  std::endl;
    for(std::size_t i = 0;  i < refine; ++i){
        ofile << i << " "<< x_top[i] << " " << y_top[i] << " " << h << std::endl;
    }
    for(std::size_t i = 0;  i < refine; ++i){
        ofile << i << " "<< x_bot[i] << " " << y_bot[i] << " " << 0.0 << std::endl;
    }

    ///write faces
    ofile << refine << " " << 0 << std::endl;
    std::size_t rm1 = refine-1;
    for (std::size_t i = 0; i < rm1; ++i){
        ofile << 4  << " " << index_top[i] << " " << index_top[i+1]  << " " << index_bot[i+1] << " " << index_bot[i] << std::endl;
    }
    ofile << 4 << " " << index_top[rm1] << " " << index_top[0] << " " << index_bot[0] << " " <<index_bot[rm1] << std::endl;
    ofile << 0 << std::endl << 0 << std::endl;
    ofile.close(); //redundant!!! destructor automatically closes ofstream
}

template<typename T, typename S>
void CylinderPLC(T radius, S height, std::size_t refine = 6){ //refine defines the number of points 
///this code is for a cylinder standing upright.
    auto pi = PI<double>;
    std::vector<double> x,y,z,theta;
    std::vector<int> index_top, index_bot;

    theta.reserve(refine); x.reserve(refine); y.reserve(refine);z.reserve(refine);
    linspace(refine, 0.0, 2*PI<double>, theta); //populate theta with values btw 0 and 2pi using linear progression

    for(std::size_t i = 0; i < refine; ++i){
        x.push_back(radius*cos(theta[i]));
        y.push_back(radius*sin(theta[i]));
        index_top.push_back(i);
        index_bot.push_back(i+refine);
    }
    write_smesh(x,y,x,y,index_top,index_bot,height,refine);
}


int main(){

    CylinderPLC(10,30,100);
    #if 0  //testing linspace
    std::vector<double> test; test.reserve(10);
    linspace(20,1,10,test);
    for(auto m:test){
        std::cout << m << std::endl;
    }
    #endif

}