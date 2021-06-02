/**
    \file poisson.cpp
    How to solve the poisson equation in a bidimensional domain using NeoPZ
*/
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGeoMeshTools.h> //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <MMeshType.h> //for MMeshType
#include <pzmanvector.h>//for TPZManVector
#include <Projection/TPZHCurlProjection.h> //for TPZHCurlProjection
#include <TPZBndCond.h> //for TPZBndCond
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>
#include <pzlog.h>

int main(int argc, char *argv[])
{
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
#ifdef PZ_LOG
   TPZLogger::InitializePZLOG();
#endif
  /* We will project an analytic solution on a HCurl approximation space
   * in the domain Omega=[-1,1]x[-1,1] embedded in a 3D space*/

  /* Suggested integration rule. since it is negative,
   * it will use the same order as the element*/
  constexpr int solOrder{-1};
  auto exactSol = [](const TPZVec<REAL> &loc,
                     TPZVec<STATE> &u,
                     TPZFMatrix<STATE> &curlU) {
    const auto &x = loc[0];
    const auto &y = loc[1];
    // u[0] = sin(M_PI * y);
    // u[1] = sin(M_PI * x);
    // curlU(0, 0) = M_PI * cos(M_PI * x) - M_PI * cos(M_PI * y);
    u[0] = (y - 1) * (y + 1);
    u[1] = (x - 1) * (x + 1);
    curlU(0, 0) = 2 * x - 2 * y;
  };
  //for the rhs, the exact sol should be sol + curl of sol
  const auto rhs = [exactSol](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
      TPZFNMatrix<1,STATE> curlU(1,1);
      exactSol(loc,u,curlU);
      u[2]=curlU(0,0);
  };
  
  //dimension of the problem
  constexpr int dim{2};
  //n divisions in x direction
  constexpr int nDivX{32};
  //n divisions in y direction
  constexpr int nDivY{32};
  
  //TPZManVector<Type,N> is a vector container with static + dynamic storage. one can also use TPZVec<Type> for dynamic storage
  TPZManVector<int,2> nDivs={nDivX,nDivY};

  //all geometric coordinates in NeoPZ are in the 3D space

  //lower left corner of the domain
  TPZManVector<REAL,3> minX={-1,-1,0};
  //upper right corner of the domain
  TPZManVector<REAL,3> maxX={ 1, 1,0};

  /*vector for storing materials(regions) identifiers
   * for the TPZGeoMeshTools::CreateGeoMeshOnGrid function.
   * In this example, we have different materials in each
   * section of the boundary*/
  TPZManVector<int,5> matIdVec={1,-1,-2,-3,-4};
  //whether to create boundary elements
  constexpr bool genBoundEls{true};

  /* val1 and val2 are used for calculating the boundary
   * conditions. val1 goes in the matrix and val2 in the rhs.
   * for dirichlet boundary conditions, only the value of
   * val2 is used.
   By default, dirichlet = 0 is imposed in all boundaries
  */
  TPZManVector<TPZFMatrix<STATE>,4> val1(4,TPZFMatrix<STATE>(1, 1, 0.));
  TPZManVector<TPZManVector<STATE, 1>,4> val2 = {{0},{0},{0},{0}};
  // dirichlet=0,neumann=1,robin=2
  TPZManVector<int,4> boundType = {0,0,0,0};
  //type of elements
  constexpr MMeshType meshType{MMeshType::ETriangular};

  //defining the geometry of the problem
  //TPZAutoPointer is a smart pointer from the NeoPZ library
  TPZAutoPointer<TPZGeoMesh> gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,genBoundEls);
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

  auto *mat = new TPZHCurlProjection<STATE>(matIdVec[0],dim);

  mat->SetForcingFunction(rhs,solOrder);
  cmesh->InsertMaterialObject(mat);

  //now we insert the boundary conditions
  for(auto i = 0; i < 4; i++)
    {
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * bnd = mat->CreateBC(mat,matIdVec[i+1],boundType[i],val1[i],val2[i]);
      cmesh->InsertMaterialObject(bnd);
    }
  
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
  //in this case, a symmetric sparse matrix is used if NeoPZ
  // was configured with MKL, otherwise, a sym. skyline matrix is used
#ifdef PZ_USING_MKL
  TPZSSpStructMatrix<STATE> strmat(cmesh);
#else
  TPZSkylineStructMatrix<STATE> strmat(cmesh);
#endif
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
    TPZSimpleTimer err("Calc error");
    // let us set the exact solution and suggest an integration rule
    // for calculating the error
    an.SetExact(exactSol, solOrder);

    /// Calculating approximation error
    
    std::ofstream anPostProcessFile("postprocess.txt");
    an.SetThreadsForError(nThreads);
    an.PostProcess(error, anPostProcessFile);
  }
	
  std::cout << "\nApproximation error:\n"
            << "HCurl Norm = " << error[0]<<'\n'
            << "L2 Norm = " << error[1]<<'\n'
            << "HCurl Seminorm = " << error[2] << "\n\n";
            
  ///vtk export
  TPZVec<std::string> scalarVars(1), vectorVars(1);
  vectorVars[0] = "Solution";
  scalarVars[0] = "Curl";
  an.DefineGraphMesh(2,scalarVars,vectorVars,"hcurlProjection.vtk");
  constexpr int resolution{1};

  std::cout << "Post procesing..."<<std::endl;
  TPZSimpleTimer post("Post-processing");
  
  an.PostProcess(resolution);	
  return 0;
}

