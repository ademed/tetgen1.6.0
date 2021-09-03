#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <utility>
#include <chrono>
#include "mesh.hpp"
#include <omp.h>
//#include "tetgen.h"




int main(){
    //mxcpl::mesh::Cylinder a(3,10,50); 

    mxcpl::mesh::Rectangular_Cuboid rec(5,10,2);
    mxcpl::mesh::tetrahedra_mesh M(rec);     
    M.output("rec");

    // #pragma omp parallel for
    // for(int i = 0; i<4; ++i)
    // {
    //     std::cout << "Ahmed" << "\n";   
    // }
    
    // std::ofstream ofile("cylinder.mtr"); 
    // ofile << 100 << " " << 1 <<"\n";
    // for(std::size_t i = 0; i < 100; ++i){
    //     ofile <<  0.125 << "\n";
    // }
    

    //test code for running tetgen library
    #if 0  
    tetgenio in, out;
    tetgenio::facet *f;
    tetgenio::polygon *p;
    int i;

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
    tetgenbehavior b; b.parse_commandline("pq1.414a0.1fnn"); b.vtkview = 1; //tetgenmesh::outmesh2vtk
    tetrahedralize(&b, &in, &out); 

    //Output mesh to files 'barout.node', 'barout.ele' and 'barout.face'.
     std::ofstream ofile("indices.txt");
    for (size_t i = 0; i < 4165; i++)
    {
        ofile << *out.face2tetlist++ << "\n"; 

    }
   std::cout <<  "faces: "<< out.numberoftrifaces << "\n" << "tetra: " << out.numberoftetrahedra << std::endl; 
    out.save_nodes("barout");
    out.save_elements("barout");
    out.save_faces("barout"); out.save_neighbors("barout"); 
    #endif
 
    #if 0
    auto refine = 5;
    //int  m = refine*(2+10); // write_smesh
    int m = refine*2;  //write_smesh2
    auto start = std::chrono::high_resolution_clock::now(); 
     CYLINDER_OBJECT c(5.753,15,refine); c.get_smesh();
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;
    std::ofstream ofile("c.mtr"); 
    ofile << m << " " << 1 <<"\n";
    for(std::size_t i = 0; i < m; ++i){
        ofile <<  0.125 << "\n";
    }
    
    
    int i = 0; /////FINALLY UNDERSTOOD THE IMPLICATIONS OF ++I AND I++
    std::vector<int> v = {1,2,3,4,5};
    for (size_t j = 0; j < 5; j++)
    {
        std::cout << v[i++] << std::endl;
    }
    #endif
}



 