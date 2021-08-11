#pragma once
#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include <fstream>
#include "mxcpl_macros.hpp"


using vector_d      =  typename std::vector<double>;
using indices       =  typename std::vector<int>; 
using u_int         =  std::size_t; 



_MXCPL_MESH_BEGIN
    template<typename T>
    constexpr T PI{ 3.14159265358979323846264338327950288419716939937510582};

    constexpr void linspace(int npoints, double first, double last, vector_d& res);

    template<typename T, typename S>
    struct Cylinder{
        public:
            Cylinder(T _radius, S _height, u_int nPoints = 10): 
                    mRadius(_radius), mHeight(_height), mNumber_of_points(nPoints){create();}
            void get_smesh(){ write_smesh();}
            

        private:
            void create();
            void write_smesh();
            T mRadius;
            S mHeight;
            u_int mNumber_of_points;
            vector_d mX, mY, mTheta; //vector of points that generate the top and bottom cirle(polygon) of cylinder
            indices mIndex_top, mIndex_bot;
    };
 

    
template<typename _T, typename _S>
void Cylinder<_T,_S>::create(){
    mTheta.reserve(mNumber_of_points); mX.reserve(mNumber_of_points); mY.reserve(mNumber_of_points);
    linspace(mNumber_of_points, 0.0, 2*PI<double> - 0.1, mTheta); //populate theta with values btw 0 and 2pi using linear progression

    for(u_int i = 0; i < mNumber_of_points; ++i){
        mX.push_back(mRadius*cos(mTheta[i]));
        mY.push_back(mRadius*sin(mTheta[i]));
        mIndex_top.push_back(i); //store indices of points of top polygon
        mIndex_bot.push_back(i+mNumber_of_points); //store indices of points of base poiygon
    }
}

template<typename _T, typename _S>
void Cylinder<_T,_S>::write_smesh(){
    std::ofstream ofile("c.smesh"); 
    //WRITE NODES-
    u_int N = 2 * mNumber_of_points;
    ofile << N << " " << 3 <<" " << 0 << " " << 0 <<  std::endl;
    for(u_int i = 0;  i < mNumber_of_points; ++i){
        ofile << mIndex_top[i] << " "<< mX[i] << " " << mY[i] << " " << mHeight << "\n";
    }
    for(u_int i = 0;  i < mNumber_of_points; ++i){
        ofile << mIndex_bot[i] << " "<< mX[i] << " " << mY[i] << " " << 0.0 << "\n";
    }

    //WRITE FACES-
    ofile << mNumber_of_points + 2 << " " << 0 << "\n" << mNumber_of_points << " ";
    //ofile << refine << " " << 0 ;

    for(u_int i = 0; i < mNumber_of_points; ++i){
        ofile << mIndex_top[i] << " "; //top polygon for cylinder
    }
    ofile << "\n" << mNumber_of_points << " ";
    for(u_int i = 0; i < mNumber_of_points; ++i){
        ofile << mIndex_bot[i] << " ";  //bottom polygon for cylinder
    }
    ofile << "\n";

    u_int rm1 = mNumber_of_points - 1;
    for (u_int i = 0; i < rm1; ++i){
        ofile << 4  << " " << mIndex_top[i] << " " << mIndex_top[i+1]  << " " << mIndex_bot[i+1] << " " << mIndex_bot[i] << "\n";
    }   
    ofile << 4 << " " << mIndex_top[rm1] << " " << mIndex_top[0] << " " << mIndex_bot[0] << " " <<mIndex_bot[rm1] << "\n"; 
    ofile << 0 << "\n" << 0 << std::endl;
    ofile.close(); //redundant!!! destructor automatically closes ofstream
}



/////////////////////////////////////////////////////////////
//                                                          //                                                    
//     Utilities                                            //                              
//                                                          //                                                    
////////////////////////////////////////////////////////////
constexpr void linspace(int npoints, double first, double last, vector_d& res){
    double dist = (last - first)/(npoints - 1); //get common difference between terms of arithmetic progression
    for (u_int i = 0; i < npoints; ++i)
    {
        res.push_back(first + i*dist);
    }   
}

_MXCPL_MESH_END

#endif //_MESH_