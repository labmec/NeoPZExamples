/**
    \file HCurlProjection.cpp
    How to project an analytic solution in a HCurl-conforming approximation space.
*/
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGmshReader.h>
#include <pzmanvector.h>//for TPZManVector
#include <TPZBndCond.h> //for TPZBndCond
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>
#include <TPZVTKGenerator.h>
#include <pzlog.h>

#include "TPZMatHCurl3D.h"//for TPZMatHCurlProjection

int main(int argc, char *argv[])
{
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
#ifdef PZ_LOG
   TPZLogger::InitializePZLOG();
#endif

   //coefficient of the mass term of the equation
   constexpr STATE coeff{1};

   //max polynomial order of the rhs
   constexpr auto rhsOrder{6};
   const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
     const auto &x = loc[0];
     const auto &y = loc[1];
     const auto &z = loc[2];
     const auto onemx2 = 1-x*x;
     const auto onemy2 = 1-y*y;
     const auto onemz2 = 1-z*z;
     u[0] = x*y*onemy2*onemz2+2*x*y*onemz2;
     u[1] = y*y*onemx2*onemz2 + onemy2*(2-x*x-z*z);
     u[2] = y*z*onemx2*onemy2 + 2*y*z*onemx2;
   };
  
   //dimension of the problem
   constexpr int dim{3};
   //whether to create boundary elements
   constexpr bool genBoundEls{true};

   /* val1 and val2 are used for calculating the boundary
    * conditions. val1 goes in the matrix and val2 in the rhs.
    * for dirichlet boundary conditions, only the value of
    * val2 is used.
    By default, dirichlet = 0 is imposed in all boundaries
   */
   TPZFMatrix<STATE> val1(1,1,0.);
   TPZManVector<STATE, 1> val2 = {0};
   // dirichlet=0,neumann=1,robin=2
   constexpr int boundType = {0};
   //type of elements

   //the geometry/mesh of this problem is defined by a .msh mesh
   const std::string filename{"hcurlmesh.msh"};
   /*this structure will give us information about the materials read by GmshReader
     position i will give us all the physical entities of dimension i as a map
     having as keys the name of the region and as values their identifier
    */
   TPZVec<std::map<std::string,int>> mat_info;
   //TPZAutoPointer is a smart pointer from the NeoPZ library
   TPZAutoPointer<TPZGeoMesh> gmesh = [&filename, &mat_info] (){
     TPZGmshReader meshReader;
     auto gmesh = meshReader.GeometricGmshMesh(filename);
     
     mat_info = meshReader.GetDimNamePhysical();
     const auto dim = meshReader.Dimension();
     for(int i = 0; i <= dim; i++){
       std::cout<<"materials with dim "<<i<<std::endl;
       for(auto &mat : mat_info[i]){
         std::cout<<"\t name: "<<mat.first <<" id: "<<mat.second<<std::endl;
       }
     }
     return gmesh;
   }();

   const int volid = mat_info[3]["vol"];
   const int bndid = mat_info[2]["bnd"];
   
   ///Defines the computational mesh based on the geometric mesh
   TPZAutoPointer<TPZCompMesh>  cmesh = new TPZCompMesh(gmesh);

   //polynomial order used in the approximatoin
   constexpr int pOrder{2};
   //using HCurl-conforming elements
   cmesh->SetAllCreateFunctionsHCurl();

   /* The TPZMaterial class is used for implementing the weak formulation.
    * Each instance has an associated material id, which should correspond to the
    * material ids used when creating the geometric mesh. In this way, you could
    * have different materials on different mesh regions */

   auto *mat = new TPZMatHCurl3D(volid,coeff);

   mat->SetForcingFunction(rhs,rhsOrder);
   cmesh->InsertMaterialObject(mat);

   //now we insert the boundary condition
   //TPZBndCond is a material type for boundary conditions
   TPZBndCond * bnd = mat->CreateBC(mat,bndid,boundType,val1, val2);
   cmesh->InsertMaterialObject(bnd);
   
   //seting the polynomial order in the computational mesh
   cmesh->SetDefaultOrder(pOrder);
   //actually creates the computational elements
   cmesh->AutoBuild();

   /*The TPZLinearAnalysis class manages the creation of the algebric
    * problem and the matrix inversion*/
   TPZLinearAnalysis an(cmesh);

   //sets number of threads to be used by the solver
   constexpr int nThreads{8};
   //defines storage scheme to be used for the FEM matrices
   //in this case, a symmetric sparse matrix is used (needs MKL)
   TPZSSpStructMatrix<STATE> strmat(cmesh);

   strmat.SetNumThreads(nThreads);
   an.SetStructuralMatrix(strmat);
  	
   ///Setting a direct solver
   TPZStepSolver<STATE> step;
   step.SetDirect(ELDLt);
   an.SetSolver(step);

   TPZManVector<REAL, 3> error;
   {
     TPZSimpleTimer total("Total");
     {
       TPZSimpleTimer assemble("Assemble");
       //assembles the system
       an.Assemble();
     }
     {
       TPZSimpleTimer solve("Solve");
       ///solves the system
       an.Solve();
     }
   }
            
   ///vtk export
   std::cout << "Post processing..."<<std::endl;
   TPZSimpleTimer tpostprocess("Post processing");
   TPZVec<std::string> fvars = {
      "u",
      "curl_u"};
   const auto file{"hcurl3d"};
   constexpr auto vtkRes{2};
   auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkRes);
   vtk.Do();
   return 0;
}

