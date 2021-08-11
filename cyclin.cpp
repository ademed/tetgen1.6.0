#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>



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
                std::vector<int> const& index_top,std::vector<int> const& index_bot, double h, std::size_t refine, std::size_t n_space = 10 ){
    
    std::ofstream ofile("c.smesh"); 
    ////----------------WRITE NODES-----------------------//////////
    std::size_t N = 2*refine;
    // ofile << N << " " << 3 <<" " << 0 << " " << 0 <<  std::endl;
    ofile << refine*(2+n_space) << " " << 3 <<" " << 0 << " " << 0 <<  std::endl;
    for(std::size_t i = 0;  i < refine; ++i){
        ofile << index_top[i] << " "<< x_top[i] << " " << y_top[i] << " " << h << "\n";
    }
    for(std::size_t i = 0;  i < refine; ++i){
        ofile << index_bot[i] << " "<< x_bot[i] << " " << y_bot[i] << " " << 0.0 << "\n";
    }
    std::vector<int> index_extra; index_extra.reserve(n_space*refine);
    double spacing = h/n_space;
    for(std::size_t i = 0; i < refine; ++i){
        for (std::size_t j = 1; j <= n_space ; ++j)
        {
            index_extra.push_back(N);
            ofile << N << " " << x_top[i]  << " " << y_top[i]  << " " << h - j*spacing << "\n";
            ++N;
        }
    }

    //////--------------WRITE FACES-------------------------------////////
    //ofile << refine + 2 << " " << 0 << "\n" << refine << " ";
     ofile << refine *(n_space + 1) + 2 << " " << 0 << "\n" << refine << " ";
    //ofile << refine << " " << 0 ;

    for(std::size_t i = 0; i < refine; ++i){
        ofile << index_top[i] << " "; //top polygon for a cylinder
    }
    ofile << "\n" << refine << " ";
    for(std::size_t i = 0; i < refine; ++i){
        ofile << index_bot[i] << " ";  //bottom polygon for a cylinder
    }
    ofile << "\n";

    std::size_t rm1 = refine-1;
    // for (std::size_t i = 0; i < rm1; ++i){
    //     ofile << 4  << " " << index_top[i] << " " << index_top[i+1]  << " " << index_bot[i+1] << " " << index_bot[i] << "\n";
    // }
    for (std::size_t i = 0; i < rm1; ++i){
        ofile << 4  << " " << index_top[i] << " " << index_top[i+1]  << " " << index_extra[(i+1)*n_space] << " " << index_extra[i*n_space] << "\n";
        for (std::size_t j = 0; j < n_space - 1; ++j) {
             ofile << 4  << " " << index_extra[j + i*n_space] << " " << index_extra[j + (i+1)*n_space] 
              << " " << index_extra[j + 1 + (i+1)*n_space] << " " << index_extra[j + 1 + i*n_space] << "\n";
        }
        ofile << 4  << " " << index_extra[(i+1)*n_space - 1] << " " << index_extra[(i+2)*n_space - 1]  << " " << index_bot[i+1] << " " << index_bot[i] << "\n";
    }
    
    ofile << 4  << " " << index_top[rm1] << " " << index_top[0]  << " " << index_extra[0] << " " << index_extra[rm1*n_space] << "\n";
        for (std::size_t j = 0; j < n_space - 1; ++j) {
             ofile << 4  << " " << index_extra[j + rm1*n_space] << " " << index_extra[j] 
              << " " << index_extra[j + 1] << " " << index_extra[j + 1 + rm1*n_space] << "\n";
        }
    ofile << 4  << " " << index_extra[n_space - 1] << " " << index_extra[refine*n_space - 1]  << " " << index_bot[0] << " " << index_bot[rm1] << "\n";
    
    //ofile << 4 << " " << index_top[rm1] << " " << index_top[0] << " " << index_bot[0] << " " <<index_bot[rm1] << "\n";
    
    
    ofile << 0 << "\n" << 0 << std::endl;
    ofile.close(); //redundant!!! destructor automatically closes ofstream
}






template<typename T>
void write_smesh2(T const& x_top, T const& y_top, T const& x_bot, T const& y_bot,
                std::vector<int> const& index_top,std::vector<int> const& index_bot, double h, std::size_t refine ){
    
    std::ofstream ofile("c.smesh"); 
    ////----------------WRITE NODES-----------------------//////////
    std::size_t N = 2*refine;
    ofile << N << " " << 3 <<" " << 0 << " " << 0 <<  std::endl;
    for(std::size_t i = 0;  i < refine; ++i){
        ofile << index_top[i] << " "<< x_top[i] << " " << y_top[i] << " " << h << "\n";
    }
    for(std::size_t i = 0;  i < refine; ++i){
        ofile << index_bot[i] << " "<< x_bot[i] << " " << y_bot[i] << " " << 0.0 << "\n";
    }

    //////--------------WRITE FACES-------------------------------////////
    ofile << refine + 2 << " " << 0 << "\n" << refine << " ";
    //ofile << refine << " " << 0 ;

    for(std::size_t i = 0; i < refine; ++i){
        ofile << index_top[i] << " "; //top polygon for a cylinder
    }
    ofile << "\n" << refine << " ";
    for(std::size_t i = 0; i < refine; ++i){
        ofile << index_bot[i] << " ";  //bottom polygon for a cylinder
    }
    ofile << "\n";

    std::size_t rm1 = refine-1;
    for (std::size_t i = 0; i < rm1; ++i){
        ofile << 4  << " " << index_top[i] << " " << index_top[i+1]  << " " << index_bot[i+1] << " " << index_bot[i] << "\n";
    }   
    ofile << 4 << " " << index_top[rm1] << " " << index_top[0] << " " << index_bot[0] << " " <<index_bot[rm1] << "\n"; 
    ofile << 0 << "\n" << 0 << std::endl;
    ofile.close(); //redundant!!! destructor automatically closes ofstream
}

template<typename T, typename S>
void CylinderPLC(T radius, S height, std::size_t refine = 10){ //refine defines the number of points 
///this code is for a cylinder standing upright.
    auto pi = PI<double>;
    std::vector<double> x,y,z,theta;
    std::vector<int> index_top, index_bot;

    theta.reserve(refine); x.reserve(refine); y.reserve(refine);z.reserve(refine);
    linspace(refine, 0.0, 2*PI<double> - 0.1, theta); //populate theta with values btw 0 and 2pi using linear progression

    for(std::size_t i = 0; i < refine; ++i){
        x.push_back(radius*cos(theta[i]));
        y.push_back(radius*sin(theta[i]));
        index_top.push_back(i); //store indices of points of top polygon
        index_bot.push_back(i+refine); //store indices of points of base poiygon
    }
    write_smesh2(x,y,x,y,index_top,index_bot,height,refine);
}


int main(){
    auto refine = 10;
   // int  m = refine*(2+10); // write_smesh
    int m = refine*2;  //write_smesh2
    auto start = std::chrono::high_resolution_clock::now();
    CylinderPLC(4,10,refine);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;
    std::ofstream ofile("c.mtr");
    ofile << m << " " << 1 <<"\n";
    for(std::size_t i = 0; i < m; ++i){
        ofile <<  0.2 << "\n";
    }

    #if 0  //testing linspace
    std::vector<double> test; test.reserve(10);
    linspace(20,1,10,test);
    for(auto m:test){
        std::cout << m << std::endl;
    }
    #endif

}