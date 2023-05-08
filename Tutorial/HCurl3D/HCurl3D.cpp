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
#include <Poisson/TPZMatPoisson.h>
#include <pzlog.h>

#include "TPZMatHCurl3D.h"//for TPZMatHCurlProjection

// #define DEBUG_POISSON

void
FilterBoundaryEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                         TPZVec<int64_t> &activeEquations,
                         std::set<int64_t> &boundConnects);

int main(int argc, char *argv[])
{
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
#ifdef PZ_LOG
   TPZLogger::InitializePZLOG();
#endif

   constexpr bool is_complex{false};

#ifdef DEBUG_POISSON
   //coefficient of the mass term of the equation
   constexpr STATE coeff{1};

   //max polynomial order of the rhs
   constexpr auto rhsOrder{6};
   const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
     const auto &x = loc[0];
     const auto &y = loc[1];
     const auto &z = loc[2];
     const auto sinx = sin(M_PI*x);
     const auto siny = sin(M_PI*y);
     const auto sinz = sin(M_PI*z);
     u[0] = M_PI*M_PI*sinx*siny*sinz;
   };
#else
   constexpr auto rhsOrder{4};
   constexpr STATE a_coeff{1};
   constexpr STATE c_coeff{1};
   
   const auto rhs = [a_coeff, c_coeff](const TPZVec<REAL>&loc, TPZVec<STATE> &f){
     const auto &x = loc[0];
     const auto &y = loc[1];
     const auto &z = loc[2];
     const STATE x2mone = x*x-1;
     const STATE y2mone = y*y-1;
     const STATE z2mone = z*z-1;
     f[0] = -2.*(a_coeff*z2mone + a_coeff*y2mone + 0.5*c_coeff * z2mone * y2mone);
     f[1] = -2.*(a_coeff*x2mone + a_coeff*z2mone + 0.5*c_coeff * x2mone * z2mone);
     f[2] = -2.*(a_coeff*y2mone + a_coeff*x2mone + 0.5*c_coeff * y2mone * x2mone);
   };

   };

  #endif
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
   TPZAutoPointer<TPZCompMesh>  cmesh = new TPZCompMesh(gmesh,is_complex);

   //polynomial order used in the approximatoin
   constexpr int pOrder{0};
   //using HCurl-conforming elements
#ifdef DEBUG_POISSON
   cmesh->SetAllCreateFunctionsContinuous();
#else
   cmesh->SetAllCreateFunctionsHCurl();
#endif
   /* The TPZMaterial class is used for implementing the weak formulation.
    * Each instance has an associated material id, which should correspond to the
    * material ids used when creating the geometric mesh. In this way, you could
    * have different materials on different mesh regions */
#ifdef DEBUG_POISSON
   auto *mat = new TPZMatPoisson<>(volid, dim);
#else
   auto *mat = new TPZMatHCurl3D<STATE>(volid,a_coeff, c_coeff);
#endif
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
   //in this case, a symmetric sparse matrix is used
   TPZSSpStructMatrix<STATE> strmat(cmesh);

   strmat.SetNumThreads(nThreads);
   const bool filter_bound{true};
   if(filter_bound){
     const int n_dofs_before = cmesh->NEquations();
     std::set<int64_t> boundConnects;
     TPZVec<int64_t> activeEquations;
     FilterBoundaryEquations(cmesh, activeEquations,boundConnects);
     const int n_dofs_after = activeEquations.size();
     std::cout<<"neq(before): "<<n_dofs_before
              <<"\tneq(after): "<<n_dofs_after<<std::endl;
     strmat.EquationFilter().SetActiveEquations(activeEquations);
   }else{
     std::cout<<"neq: "<<cmesh->NEquations()<<std::endl;;
   }
   
   an.SetStructuralMatrix(strmat);
  

   TPZStepSolver<STATE> solver;
   solver.SetDirect(ELDLt);
   an.SetSolver(solver);
   
   {
     TPZSimpleTimer total("Total");
     {
       TPZSimpleTimer assemble("Assemble");
       //assembles the system
       an.Assemble();
     }
     {
       TPZSimpleTimer solve("Solve", true);

       auto *solver = dynamic_cast<TPZStepSolver<STATE> *>(an.Solver());
       solver->Matrix()->SetSymmetry(SymProp::Sym);

       Precond::Type precond_type = Precond::Element;
       constexpr bool overlap {false};
       TPZAutoPointer<TPZMatrixSolver<STATE>> precond =
         an.BuildPreconditioner<STATE>(precond_type, overlap);

       
       

       const int64_t n_iter = {500};
       const int n_vecs = {150};
       constexpr REAL tol = 1e-10;
       constexpr int64_t from_current{0};
       solver->SetGMRES(n_iter, n_vecs, *precond, tol, from_current);
       //solves the system
       an.Solve();
     }
   }
   
   ///vtk export
   std::cout << "Post processing..."<<std::endl;
   TPZSimpleTimer tpostprocess("Post processing");
#ifdef DEBUG_POISSON
   TPZVec<std::string> fvars = {
      "Solution",
      "Derivative"};
#else
   TPZVec<std::string> fvars = {
      "u",
      "curl_u"};
#endif
   const auto file{"hcurl3d"};
   constexpr auto vtkRes{2};
   auto vtk = TPZVTKGenerator(cmesh, fvars, file, vtkRes);
   vtk.Do();
   return 0;
}

void
FilterBoundaryEquations(TPZAutoPointer<TPZCompMesh> cmesh,
                         TPZVec<int64_t> &activeEquations,
                         std::set<int64_t> &boundConnects)
{
  TPZSimpleTimer timer ("Filter dirichlet eqs");
  TPZManVector<int64_t, 1000> allConnects;
  boundConnects.clear();

  for (int iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if (cel == nullptr) {
      continue;
    }
    if (cel->Reference() == nullptr) {//there is no associated geometric el
      continue;
    }
    TPZBndCond *mat = dynamic_cast<TPZBndCond *>(
        cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
    if (mat && mat->Type() == 0) {//check for dirichlet bcs
      std::set<int64_t> boundConnectsEl;
      cel->BuildConnectList(boundConnectsEl);
      for(auto val : boundConnectsEl){
        if (boundConnects.find(val) == boundConnects.end()) {
          boundConnects.insert(val);
        }
      }
    }
  }

  //certainly we have less equations than this, but we will avoid repeated resizes
  activeEquations.Resize(cmesh->NEquations());
  int neq = 0;
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (boundConnects.find(iCon) == boundConnects.end()) {
      TPZConnect &con = cmesh->ConnectVec()[iCon];
      const auto hasdep = con.HasDependency();
      const auto seqnum = con.SequenceNumber();
      const auto pos = cmesh->Block().Position(seqnum);
      const auto blocksize = cmesh->Block().Size(seqnum);
      
      if(hasdep || seqnum < 0 || !blocksize) { continue; }
      const auto vs = neq;
      for (auto ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
      neq += blocksize;
    }
  }
  activeEquations.Resize(neq);
}