#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <chrono>
#include "mesh.hpp"
<<<<<<< HEAD
#include "tetgen.h"
=======
>>>>>>> 464193039a1791db0d7e708b0d14f27a6226f5bd




int main(int argc, char *argv[]){
        
    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *p;
    int i;

<<<<<<< HEAD
    // All indices start from 1.
    in.firstnumber = 1;

    in.numberofpoints = 8;
    in.pointlist = new REAL[in.numberofpoints * 3];
    in.pointlist[0]  = 0;  // node 1.
    in.pointlist[1]  = 0;
    in.pointlist[2]  = 0;
    in.pointlist[3]  = 2;  // node 2.
    in.pointlist[4]  = 0;
    in.pointlist[5]  = 0;
    in.pointlist[6]  = 2;  // node 3.
    in.pointlist[7]  = 2;
    in.pointlist[8]  = 0;
    in.pointlist[9]  = 0;  // node 4.
    in.pointlist[10] = 2;
    in.pointlist[11] = 0;
    // Set node 5, 6, 7, 8.
    for (i = 4; i < 8; i++) {
        in.pointlist[i * 3]     = in.pointlist[(i - 4) * 3];
        in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
        in.pointlist[i * 3 + 2] = 12;
    }

    in.numberoffacets = 6;
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    // Facet 1. The leftmost facet.
    f = &in.facetlist[0];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 1;
    p->vertexlist[1] = 2;
    p->vertexlist[2] = 3;
    p->vertexlist[3] = 4;
    
    // Facet 2. The rightmost facet.
    f = &in.facetlist[1];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 5;
    p->vertexlist[1] = 6;
    p->vertexlist[2] = 7;
    p->vertexlist[3] = 8;

    // Facet 3. The bottom facet.
    f = &in.facetlist[2];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 1;
    p->vertexlist[1] = 5;
    p->vertexlist[2] = 6;
    p->vertexlist[3] = 2;

    // Facet 4. The back facet.
    f = &in.facetlist[3];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 2;
    p->vertexlist[1] = 6;
    p->vertexlist[2] = 7;
    p->vertexlist[3] = 3;

    // Facet 5. The top facet.
    f = &in.facetlist[4];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 3;
    p->vertexlist[1] = 7;
    p->vertexlist[2] = 8;
    p->vertexlist[3] = 4;

    // Facet 6. The front facet.
    f = &in.facetlist[5];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0] = 4;
    p->vertexlist[1] = 8;
    p->vertexlist[2] = 5;
    p->vertexlist[3] = 1;

    // Set 'in.facetmarkerlist'

    in.facetmarkerlist[0] = -1; 
    in.facetmarkerlist[1] = -2;
    in.facetmarkerlist[2] = -3;
    in.facetmarkerlist[3] = -3;
    in.facetmarkerlist[4] = -3;
    in.facetmarkerlist[5] = -3;

    // Output the PLC to files 'barin.node' and 'barin.poly'.
    in.save_nodes("barin");
    in.save_poly("barin");

    // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
    //   do quality mesh generation (q) with a specified quality bound
    //   (1.414), and apply a maximum volume constraint (a0.1).
    tetgenbehavior b; b.parse_commandline("Qpq1.414a0.1f"); b.vtkview = 1; //tetgenmesh::outmesh2vtk
    tetrahedralize(&b, &in, &out); 

    // Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
    //  std::ofstream ofile("indices.txt");
    // for (size_t i = 0; i < 4165; i++)
    // {
    //     ofile << *out.face2tetlist++ << "\n"; 

    // }
   // std::cout <<  "faces: "<< out.numberoftrifaces << "\n" << "tetra: " << out.numberoftetrahedra << std::endl; 
    out.save_nodes("barout");
    out.save_elements("barout");
    out.save_faces("barout"); out.save_neighbors("barout"); 
 
    




  





    // auto refine = 5;
    // //int  m = refine*(2+10); // write_smesh
    // int m = refine*2;  //write_smesh2
    // auto start = std::chrono::high_resolution_clock::now(); 
    //  CYLINDER_OBJECT c(5.753,15,refine); c.get_smesh();
    // auto end = std::chrono::high_resolution_clock::now();
    // std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;
    // std::ofstream ofile("c.mtr"); 
    // ofile << m << " " << 1 <<"\n";
    // for(std::size_t i = 0; i < m; ++i){
    //     ofile <<  0.4 << "\n";
    // }
    
    
    // int i = 0; /////FINALLY UNDERSTOOD THE IMPLICATIONS OF ++I AND I++
    // std::vector<int> v = {1,2,3,4,5};
    // for (size_t j = 0; j < 5; j++)
    // {
    //     std::cout << v[i++] << std::endl;
    // }




    
=======
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
        ofile << index_top[i] << " "; //top polygon for cylinder
    }
    ofile << "\n" << refine << " ";
    for(std::size_t i = 0; i < refine; ++i){
        ofile << index_bot[i] << " ";  //bottom polygon for cylinder
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
    write_smesh(x,y,x,y,index_top,index_bot,height,refine);
}


int main(){
    auto refine = 100;
    //int  m = refine*(2+10); // write_smesh
    int m = refine*2;  //write_smesh2
    auto start = std::chrono::high_resolution_clock::now();
    CYLINDER_OBJECT c(5.753,15,refine); c.get_smesh();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;
    std::ofstream ofile("c.mtr");
    ofile << m << " " << 1 <<"\n";
    for(std::size_t i = 0; i < m; ++i){
        ofile <<  0.4 << "\n";
    }

>>>>>>> 464193039a1791db0d7e708b0d14f27a6226f5bd
    #if 0  //testing linspace
    std::vector<double> test; test.reserve(10);
    
    linspace(20,1,10,test);
    for(auto m:test){
        std::cout << m << std::endl;
    }
    
    #endif

}