/**
    \file HCurlProjection.cpp
    How to project an analytic solution in a HCurl-conforming approximation space.
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

#include "TPZShapeHDivKernel.h"
#include "pzshapetetra.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzshapecube.h"
#include "pzshapeprism.h"
#include "TPZShapeData.h"
#include "pzquad.h"


template<class TSHAPE>
void FindHCurlDependency(int order);

typedef std::pair<MElementType,int> orderpair;
std::map<orderpair ,std::set<int>> ShapeRemove;


int main(int argc,char *argv[])
{
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
#ifdef PZ_LOG
   TPZLogger::InitializePZLOG();
#endif
    ShapeRemove[orderpair(ETriangle,2)] = {0};
    ShapeRemove[orderpair(ETriangle,3)] = {0,1,7};
    ShapeRemove[orderpair(ETriangle,4)] = {6,9,10,11,13,14};
    ShapeRemove[orderpair(ETriangle,5)] = {8,12,13,14,16,17,20,21,22,23};
    ShapeRemove[orderpair(ETetraedro,3)] = {0};
    ShapeRemove[orderpair(ETetraedro,4)] = {9,12,13,14};
    ShapeRemove[orderpair(ETetraedro,5)] = {18,24,25,26,29,31,32,33,34,35};
    ShapeRemove[orderpair(EQuadrilateral,1)] = {3};
    ShapeRemove[orderpair(EQuadrilateral,2)] = {6,8,9,11};
    ShapeRemove[orderpair(EQuadrilateral,3)] = {9,12,13,15,17,18,20,22,23};
    ShapeRemove[orderpair(EQuadrilateral,4)] = {12,16,17,19,21,23,24,26,28,30,31,33,35,37,38,39};
    ShapeRemove[orderpair(EQuadrilateral,5)] = {15,20,21,23,25,27,29,30,32,34,36,38,39,41,43,45,47,48,50,52,54,56,57,58,59};
    ShapeRemove[orderpair(ECube,1)] = {5};
    ShapeRemove[orderpair(ECube,2)] = {11,17,19,20,21,24,25,26,27,28,29,30,31,32,33,34,35};
    ShapeRemove[orderpair(ECube,3)] = {22,23,24,25,26,37,38,40,41,42,43,44,45,46,47,50,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,
                                       71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107};
    ShapeRemove[orderpair(ECube,4)] = {19,23,27,31,37,38,40,41,42,43,46,47,51,55,59,63,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
                                       86,87,90,91,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,
                                       120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,
                                       148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
                                       176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,
                                       204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239};                                    
    ShapeRemove[orderpair(EPrisma,1)] = {};
    ShapeRemove[orderpair(EPrisma,2)] = {7,8};
    ShapeRemove[orderpair(EPrisma,3)] = {21,24,25,26,29,32,33,34,35};
    ShapeRemove[orderpair(EPrisma,4)] = {5,38,39,40,41,42,44,45,46,47,48,49,50,52,53,55,56,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,
                                         75,76,77,78,79,80,81,82,83,84,85,86,87,88,89};
    ShapeRemove[orderpair(EPrisma,5)] = {5,6,8,9,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,
                                         86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,
                                         116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
                                         144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179};
    // FindHCurlDependency<pzshape::TPZShapeTriang>(5);
    FindHCurlDependency<pzshape::TPZShapeTetra>(3);
    // FindHCurlDependency<pzshape::TPZShapeQuad>(3);
    // FindHCurlDependency<pzshape::TPZShapeCube>(4);
    // FindHCurlDependency<pzshape::TPZShapePrism>(2);
    return 0;
  /* We will project an analytic solution on a HCurl approximation space
   * in the domain Omega=[-1,1]x[-1,1] embedded in a 3D space*/
  constexpr int solOrder{2};
  auto exactSol = [](const TPZVec<REAL> &loc,
                     TPZVec<STATE> &u,
                     TPZFMatrix<STATE> &curlU) {
    const auto &x = loc[0];
    const auto &y = loc[1];
    u[0] = sin(M_PI * y);
    u[1] = sin(M_PI * x);
    curlU(0,0) = M_PI * cos(M_PI * x) - M_PI * cos(M_PI * y);
    // u[0] = (y - 1) * (y + 1);
    // u[1] = (x - 1) * (x + 1);
    // curlU(0,0) = 2 * x - 2 * y;
  };
  //for the rhs,the exact sol should be sol + curl of sol
  const auto rhs = [exactSol](const TPZVec<REAL>&loc,TPZVec<STATE> &u){
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
  TPZManVector<REAL,3> maxX={ 1,1,0};

  /*vector for storing materials(regions) identifiers
   * for the TPZGeoMeshTools::CreateGeoMeshOnGrid function.
   * In this example,we have different materials in each
   * section of the boundary*/
  TPZManVector<int,5> matIdVec={1,-1,-2,-3,-4};
  //whether to create boundary elements
  constexpr bool genBoundEls{true};

  /* val1 and val2 are used for calculating the boundary
   * conditions. val1 goes in the matrix and val2 in the rhs.
   * for dirichlet boundary conditions,only the value of
   * val2 is used.
   By default,dirichlet = 0 is imposed in all boundaries
  */
  TPZManVector<TPZFMatrix<STATE>,4> val1(4,TPZFMatrix<STATE>(1,1,0.));
  TPZManVector<TPZManVector<STATE,1>,4> val2 = {{0},{0},{0},{0}};
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
   * Each instance has an associated material id,which should correspond to the
   * material ids used when creating the geometric mesh. In this way,you could
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
  //in this case,a symmetric sparse matrix is used if NeoPZ
  // was configured with MKL,otherwise,a sym. skyline matrix is used
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

  TPZManVector<REAL,3> error;
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
    an.SetExact(exactSol,solOrder);

    /// Calculating approximation error
    
    std::ofstream anPostProcessFile("postprocess.txt");
    an.SetThreadsForError(nThreads);
    an.PostProcess(error,anPostProcessFile);
  }
	
  std::cout << "\nApproximation error:\n"
            << "HCurl Norm = " << error[0]<<'\n'
            << "L2 Norm = " << error[1]<<'\n'
            << "HCurl Seminorm = " << error[2] << "\n\n";
            
  ///vtk export
  TPZVec<std::string> scalarVars(1),vectorVars(1);
  vectorVars[0] = "Solution";
  scalarVars[0] = "Curl";
  an.DefineGraphMesh(2,scalarVars,vectorVars,"hcurlProjection.vtk");
  constexpr int resolution{1};

  std::cout << "Post processing..."<<std::endl;
  TPZSimpleTimer post("Post-processing");
  
  an.PostProcess(resolution);	
  return 0;
}


template<class TSHAPE>
void FindHCurlDependency(int order)
{
    constexpr int NNodes = TSHAPE::NCornerNodes;
    constexpr int NSides = TSHAPE::NSides;
    constexpr int NConnects = NSides-NNodes;
    constexpr int dim = TSHAPE::Dimension;
    TPZVec<int64_t> nodeids(NNodes);
    for(int i=0; i<NNodes; i++) nodeids[i] = i;
    TPZVec<int> orders(NConnects,order);
    TPZShapeData data;
    TPZShapeHDivKernel<TSHAPE>::Initialize(nodeids,orders,data);
    auto nshape = TPZShapeHDivKernel<TSHAPE>::NHCurlShapeF(data);
    TPZFMatrix<REAL> phi(dim,nshape),curlphi(3,nshape);
    if(dim == 2) curlphi.Resize(1,nshape);
    TPZManVector<REAL,3> pt(dim,0.3);
    TPZShapeHDivKernel<TSHAPE>::Shape(pt,data,phi,curlphi);
    // int lastconnectsize = TPZShapeHDivKernel<TSHAPE>::NConnectShapeF(NConnects-1,data);
    int lastconnectsize = TPZShapeHDivKernel<TSHAPE>::ComputeNConnectShapeF(NConnects-1,order);
    int nH1space = TPZShapeHDivKernel<TSHAPE>::NH1ShapeF(data);
    TPZManVector<std::set<int>,20 > shapetovec(nH1space);
    for(auto it : data.fVecShapeIndex)
    {
        int shape = it.second;
        int vec = it.first;
        shapetovec[it.second].insert(it.first);
    }
    for (int ish=0; ish<nH1space; ish++) {
        std::cout << "shape func " << ish << std::endl;
        for (auto ivec : shapetovec[ish]) {
            std::cout << "vec = " << ivec << " dir = ";
            for(int i=0; i<dim; i++) std::cout << data.fMasterDirections(i,ivec) << " ";
            std::cout << std::endl;
        }
    }
    
    typename TSHAPE::IntruleType intrule(2*order);
    int npoints = intrule.NPoints();
    TPZFMatrix<REAL> L2(lastconnectsize,lastconnectsize,0.),CC(lastconnectsize,lastconnectsize,0.),
        CCR(lastconnectsize,lastconnectsize,0.);
    std::set<int> remove;
    if(ShapeRemove.find(orderpair(TSHAPE::Type(),order)) != ShapeRemove.end())
    {
        remove = ShapeRemove[orderpair(TSHAPE::Type(),order)];
        CCR.Resize(lastconnectsize-remove.size(),lastconnectsize-remove.size());
    }
    REAL weight;
    for (int ip = 0; ip<npoints; ip++) {
        intrule.Point(ip,pt,weight);
        TPZShapeHDivKernel<TSHAPE>::Shape(pt,data,phi,curlphi);

        int icount = 0;
        for (int ish=0; ish<lastconnectsize; ish++) {
            int i = nshape-lastconnectsize+ish;
            int inotremoved = true;
            if(remove.find(ish) != remove.end()) inotremoved = false;
            int jcount = 0;
            for(int jsh=0; jsh<lastconnectsize; jsh++)
            {
                int jnotremoved = true;
                if(remove.find(jsh) != remove.end()) jnotremoved = false;
                int j = nshape-lastconnectsize+jsh;
                for (int d=0; d<dim; d++) {
                    L2(ish,jsh) += phi(d,i)*phi(d,j)*weight;
                }
                if(dim == 3)
                {
                    for (int d=0; d<3; d++) {
                        CC(ish,jsh) += curlphi(d,i)*curlphi(d,j)*weight;
                        if(inotremoved && jnotremoved) CCR(icount,jcount) += curlphi(d,i)*curlphi(d,j)*weight;
                    }
                } else{
                    CC(ish,jsh) += curlphi(0,ish)*curlphi(0,jsh)*weight;
                    if(inotremoved && jnotremoved) CCR(icount,jcount) += curlphi(0,i)*curlphi(0,j)*weight;
                }
                if(jnotremoved) jcount++;
            }
            if(inotremoved) icount++;
        }
    }
    std::ofstream out("CurlCheck.nb");
    L2.Print("L2Mat = ",out,EMathematicaInput);
    CC.Print("CurlCurlMat = ",out,EMathematicaInput);
    CCR.Print("CurlSelectMat = ",out,EMathematicaInput);
}

