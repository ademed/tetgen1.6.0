
#ifndef MESH_HPP
#define MESH_HPP
#include <vector>
#include <fstream>
#include <type_traits>
#include <iostream>
#include "mxcpl_macros.hpp"
#include "tetgen.h"


using vector_d      =  typename std::vector<double>;
using indices       =  typename std::vector<int>; 
using u_int         =  std::size_t; 



_MXCPL_MESH_BEGIN
    template<typename T>
    constexpr T PI{ 3.14159265358979323846264338327950288419716939937510582};
    enum OBJECTS{CYLINDER = 1, RECTANGLE = 2} obj;

    template<typename obj> //enable if cylinder or rectangle object. I only implemented type() member funciton for these two objects
    using Enable_If_CR  = std::enable_if_t< std::is_same_v<std::void_t<decltype(std::declval<obj>().type())>,void> >;

    constexpr void linspace(int npoints, double first, double last, vector_d& res);

    struct Cylinder{
        public:
            Cylinder(double _radius, double _height, u_int nPoints = 5): 
                    mRadius(_radius), mHeight(_height), mNumber_of_points(nPoints){create();}
            Cylinder(Cylinder const& rhs): 
                mRadius(rhs.mRadius), mHeight(rhs.mHeight), mNumber_of_points(rhs.mNumber_of_points){create();}
            void get_smesh(const char* filename){ write_smesh(filename);}
            auto type(){return CYLINDER;}
            vector_d::pointer x_cord() {return mX.data();}
            vector_d::pointer y_cord() {return mY.data();}
            u_int numPoints() {return mNumber_of_points;}
            double getHeight() {return mHeight;}
            indices::pointer index_top() {return mIndex_top.data();}
            indices::pointer index_bot() {return mIndex_bot.data();}
               
        private:
            void create();
            void write_smesh(const char* filename);
            double mRadius;
            double mHeight;
            u_int mNumber_of_points;
            vector_d mX, mY, mTheta; //vector of points that generate the top and bottom cirle(polygon) of cylinder
            indices mIndex_top, mIndex_bot;
    };

    struct Rectangular_Cuboid{
        public:
            Rectangular_Cuboid(double _x, double _y, double _z):
                         mWidth(_x), mLength(_y), mHeight(_z){create();}
            auto type(){return RECTANGLE;}
            vector_d::pointer x_cord() {return mX.data();}
            vector_d::pointer y_cord() {return mY.data();}
            u_int numPoints() {return mNumber_of_points;}
            double getHeight() {return mHeight;}
            indices::pointer index_top() {return mIndex_top.data();}
            indices::pointer index_bot() {return mIndex_bot.data();}
        private:
            void create();
            double mWidth,mLength,mHeight;
            u_int mNumber_of_points;
            vector_d mX, mY; //vector of points that generate the top and bottom rectangle
            indices mIndex_top, mIndex_bot;
    };


    struct tetrahedra_mesh{
        public:
            template<typename object, typename = Enable_If_CR<object>>
            tetrahedra_mesh(object&& obj);

            void output(const char* filename);

        private:
            tetgenio in, out;
            tetgenio::facet *f;
            tetgenio::polygon *p;            
    };


/////////////////////////////////////////////////////////////
//                                                         //                                                    
//     Tetrahedra_Mesh Implementation                      //                                                          
//                                                         //                                                    
//////////////////////////////////////////////////////////////
 template<typename object, typename>
tetrahedra_mesh::tetrahedra_mesh(object&& obj){
    in.firstnumber = 0;
    u_int nPoints = obj.numPoints();
    in.numberofpoints = 2*nPoints; u_int k = 0;
    in.pointlist = new REAL[in.numberofpoints * 3];
    auto ptr_x = obj.x_cord();
    auto ptr_y = obj.y_cord();
    for (u_int i = 0; i < nPoints; i++)
    {
        in.pointlist[k] = *ptr_x;
        in.pointlist[k +  3*nPoints] = *ptr_x++;
        in.pointlist[++k] = *ptr_y;
        in.pointlist[k +  3*nPoints] = *ptr_y++;
        in.pointlist[++k] = obj.getHeight();
        in.pointlist[k + 3*nPoints] = 0.0; k++;
    }

    in.numberoffacets = nPoints + 2;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets]; //no need for markers
    for (u_int i = 0; i < in.numberoffacets; ++i){
         in.facetmarkerlist[i] = -1;
    }

    //top facet of cylinder/rectangle
    f = &in.facetlist[0];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = nPoints;
    p->vertexlist = new int[p->numberofvertices];
    auto ptr_idx_top = obj.index_top();
    for(u_int i = 0; i < nPoints; ++i){
        p->vertexlist[i] = *ptr_idx_top++;
    } 

    //bottom facet of cylinder/rectangle
    f = &in.facetlist[1];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = nPoints;
    p->vertexlist = new int[p->numberofvertices];
    auto ptr_idx_bot = obj.index_bot();
    for(u_int i = 0; i < nPoints; ++i){
        p->vertexlist[i] = *ptr_idx_bot++;
    }

    //rest of the facets but last
    u_int rm1 = nPoints - 1;
    for (u_int i = 0; i <  rm1; i++)
    {
        f = &in.facetlist[i + 2];
        f->numberofpolygons = 1;
        f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
        f->numberofholes = 0;
        f->holelist = NULL;
        p = &f->polygonlist[0];
        p->numberofvertices = 4;
        p->vertexlist = new int[p->numberofvertices];
        p->vertexlist[0] = obj.index_top()[i];
        p->vertexlist[1] = obj.index_top()[i + 1];
        p->vertexlist[2] = obj.index_bot()[i + 1];
        p->vertexlist[3] = obj.index_bot()[i];
    }
    
    //last face
    f = &in.facetlist[in.numberoffacets - 1];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = obj.index_top()[rm1];
    p->vertexlist[1] = obj.index_top()[0];
    p->vertexlist[2] = obj.index_bot()[0];
    p->vertexlist[3] = obj.index_bot()[rm1];

}

void tetrahedra_mesh::output(const char* filename){
    in.save_nodes(filename); 
    in.save_poly(filename); 
    tetgenbehavior b; b.parse_commandline("pq1.404a0.2fnn"); b.vtkview = 1; 
    tetrahedralize(&b, &in, &out); 
    out.save_nodes(filename); 
    out.save_elements(filename);
    out.save_faces(filename);
    std::ofstream ofile("connected_cells.txt");
   
    for (size_t i = 0; i <  out.numberoftrifaces; i++)
    {
        ofile << *out.face2tetlist++ << " " << *out.face2tetlist++ << "\n"; 
    }

   std::cout <<  "faces: "<< out.numberoftrifaces << "\n" << "tetra: " << out.numberoftetrahedra << std::endl; 
}

 /////////////////////////////////////////////////////////////
//                                                         //                                                    
//   Cylinder Struct Implementation                        //                                               
//                                                         //                                                    
//////////////////////////////////////////////////////////////   

void Cylinder::create(){
    mTheta.reserve(mNumber_of_points); mX.reserve(mNumber_of_points); mY.reserve(mNumber_of_points);
    linspace(mNumber_of_points, 0.0, 2*PI<double> - 0.1, mTheta); //populate theta with values btw 0 and 2pi using linear progression

    for(u_int i = 0; i < mNumber_of_points; ++i){
        mX.push_back(mRadius*cos(mTheta[i]));
        mY.push_back(mRadius*sin(mTheta[i]));
        mIndex_top.push_back(i); //store indices of points of top polygon
        mIndex_bot.push_back(i+mNumber_of_points); //store indices of points of base poiygon
    }
}

void Cylinder::write_smesh(const char* filename){
    std::ofstream ofile(filename); 
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
//                                                         //                                                    
//   Rectangular cuboid Struct Implementation              //                                               
//                                                         //                                                    
////////////////////////////////////////////////////////////// 

void Rectangular_Cuboid::create(){
    mNumber_of_points = 4;
    mX.reserve(mNumber_of_points); mY.reserve(mNumber_of_points);
    mX.push_back(mWidth);  mX.push_back(mWidth); mX.push_back(0); mX.push_back(0); //counter clockwise
    mY.push_back(0);  mY.push_back(mLength); mY.push_back(mLength); mY.push_back(0); //counter clockwise

    for(u_int i = 0; i < mNumber_of_points; ++i){
        mIndex_top.push_back(i); //store indices of points of top rectangle
        mIndex_bot.push_back(i+mNumber_of_points); //store indices of points of base rectangle
    }
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