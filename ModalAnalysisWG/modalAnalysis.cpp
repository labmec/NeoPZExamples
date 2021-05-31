#include <Electromagnetics/TPZWaveguideModalAnalysis.h> //for TPZMatWaveguideModalAnalysis
#include <MMeshType.h>                                  //for MMeshType
#include <TPZBndCond.h>                                 //for TPZBndCond
#include <TPZEigenAnalysis.h>                           //for TPZLinearAnalysis
#include <TPZGeoMeshTools.h>      //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <TPZKrylovEigenSolver.h> //for TPZKrylovEigenSolver
#include <TPZSpectralTransform.h> //for TPZSTShiftInvert
#include <TPZNullMaterial.h>      //for TPZNullMaterial
#include <TPZElectromagneticConstants.h>
#include <TPZSimpleTimer.h>
#include <TPZSkylineNSymStructMatrix.h> //non-symmetric skyline matrix storage
#include <TPZSpStructMatrix.h>          //non-symmetric sparse matrix storage
#include <pzcmesh.h>                    //for TPZCompMesh
#include <pzgmesh.h>                    //for TPZGeoMesh
#include <pzmanvector.h>                //for TPZManVector
#include "TPZVTKGeoMesh.h" //for exporting geomesh to vtk
#include <pzbuildmultiphysicsmesh.h>
#include <pzlog.h>
//! BC to be used if splitting domain in half (width-wise)
enum class ESymType { NONE, PMC, PEC };


/**
   @brief Creates the geometrical mesh associated with a rectangular metallic waveguide. 
The domain can be divided in half (width-wise) for exploring symmetry. */
TPZAutoPointer<TPZGeoMesh>
CreateGMeshRectangularWaveguide(TPZVec<int> &matIdVec, const REAL &scale,
                                const REAL wDomain, const REAL hDomain,
                                TPZVec<int> &nDivs,
                                const MMeshType meshType,
                                const bool usingSymmetry);

TPZVec<TPZAutoPointer<TPZCompMesh>>
CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
            const TPZVec<int> &matIdVec, const CSTATE ur, const CSTATE er,
            const STATE lambda, const REAL &scale, bool usingSymmetry, ESymType sym);

void FilterBoundaryEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal);

int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif
  // width of the waveguide
  constexpr REAL wDomain{0.02286};
  // height of the waveguide
  constexpr REAL hDomain{0.01016};
  //whether to split the domain in two for exploring symmetries
  constexpr bool usingSymmetry{false};
  //which BC to apply if symmetry is used
  constexpr ESymType sym{ESymType::PEC};
  // n divisions in x direction
  constexpr int nDivX{40};
  // n divisions in y direction
  constexpr int nDivY{20};
  // type of elements
  constexpr MMeshType meshType{MMeshType::ETriangular};
  // TPZManVector<Type,N> is a vector container with static + dynamic storage.
  // one can also use TPZVec<Type> for dynamic storage
  TPZManVector<int, 2> nDivs = {nDivX, nDivY};
  // polynomial order to be used in the approximation
  constexpr int pOrder{2};
  // magnetic permeability
  constexpr CSTATE ur{1};
  // electric permittivity
  constexpr CSTATE er{1};
  // operational frequency
  constexpr STATE freq = 25e9;
  // wavelength
  constexpr STATE lambda{pzeletromag::cZero/freq};
  /*Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Uing k0 as a scale factor results in 
    the eigenvalues propagationConstant/k0 = effectiveIndex*/
  constexpr REAL scale{2*M_PI/lambda};
  //number of threads to use
  constexpr int nThreads{4};
  //number of genvalues to be computed
  constexpr int nEigenpairs{10};
  //whether to compute eigenvectors
  constexpr bool computeVectors{true};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff)*/
  constexpr CSTATE target = -0.93120625;
  // Dimension of the krylov space to be used. Suggested to be at least nev * 10
  constexpr int krylovDim{200};
  
  /* Vector for storing materials(regions) identifiers for the 
     TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     In this example, we have one material for the interior of the waveguide
     and one for each BC*/
  TPZManVector<int, 5> matIdVec({1,-1,-2,-3,-4});

  TPZSimpleTimer total("Total");
  //creates geometric mesh
  auto gmesh = CreateGMeshRectangularWaveguide(matIdVec,scale,wDomain,hDomain,nDivs,
                                               meshType,usingSymmetry);

  const std::string gmeshFileName{"gmesh_" +
    std::to_string(nDivX) + " _" + std::to_string(nDivY) + ".vtk"};
  //print gmesh to vtk
  {
    std::ofstream gmeshFile(gmeshFileName);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, gmeshFile, true);
  }

  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = CreateCMesh(gmesh,pOrder,matIdVec,ur,er,lambda,
                             scale,usingSymmetry,sym);
  //geths the multiphysics mesh (main mesh)
  auto cmesh = meshVec[0];

  //reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth{true};
  TPZEigenAnalysis an(cmesh, optimizeBandwidth);
  an.SetComputeEigenvectors(computeVectors);

  /**
     When using NeoPZ with MKL, an sparse matrix should be used for better
     performance. Otherwise, the skyline matrix has available inhouse solvers.
  */
  TPZAutoPointer<TPZStructMatrix> strmtrx{nullptr};
#ifdef PZ_USING_MKL
  strmtrx = new TPZSpStructMatrix<CSTATE>(cmesh);
#else
  strmtrx = new TPZSkylineNSymStructMatrix<CSTATE>(cmesh);
#endif
  
  strmtrx->SetNumThreads(nThreads);
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system
   */
  TPZVec<int64_t> activeEquations;
  //this value is the total number of dofs including dirichlet bcs
  int neqOriginal;
  {
    int neq;
    FilterBoundaryEquations(meshVec, activeEquations, neq, neqOriginal);
  }
  strmtrx->EquationFilter().SetActiveEquations(activeEquations);
  an.SetStructuralMatrix(strmtrx);

  TPZSTShiftAndInvert<CSTATE> st;
  st.SetShift(target);
  TPZKrylovEigenSolver<CSTATE> solver;
  solver.SetSpectralTransform(st);
  solver.SetKrylovDim(krylovDim);
  solver.SetNEigenpairs(nEigenpairs);
  solver.SetAsGeneralised(true);
  
  an.SetSolver(solver);
  {
    TPZSimpleTimer assemble("Assemble");
    an.Assemble();
  }
  an.Solve();
  auto ev = an.GetEigenvalues();

  for(auto &w : ev){
    std::cout<<w<<std::endl;
  }
  
  if (!computeVectors) return 0;
  
  TPZStack<std::string> scalnames, vecnames;
  scalnames.Push("Ez");
  vecnames.Push("Et");
  const std::string plotfile = "fieldPlot.vtk";
  constexpr int dim{2};
  an.DefineGraphMesh(dim, scalnames, vecnames,plotfile);

  auto eigenvectors = an.GetEigenvectors();

  TPZFMatrix<CSTATE> evector(neqOriginal, 1, 0.);
  const auto nev = ev.size();

  TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
  meshVecPost[0] = meshVec[1];
  meshVecPost[1] = meshVec[2];
  constexpr int vtkRes{1};
  std::cout<<"Post processing..."<<std::endl;
  for (int iSol = 0; iSol < ev.size(); iSol++) {
    const CSTATE currentKz = std::sqrt(-1. * ev[iSol]);
    eigenvectors.GetSub(0, iSol, neqOriginal, 1, evector);
    for(auto id : matIdVec){
      auto matPtr =
        dynamic_cast<TPZWaveguideModalAnalysis *>(cmesh->FindMaterial(id));
      if(!matPtr) continue;
      matPtr->SetKz(currentKz);
      matPtr->SetPrintFieldRealPart(true);
    }
    std::cout<<"\rPost processing step "<<iSol<<" out of "<<ev.size()<<std::flush;
    an.LoadSolution(evector);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, cmesh);
    an.PostProcess(vtkRes);
  }
  std::cout<<"\rFinished post processing"<<std::endl;
  std::cout<<std::endl;
  return 0;
}

TPZAutoPointer<TPZGeoMesh> CreateGMeshRectangularWaveguide(
    TPZVec<int> &matIdVec, const REAL &scale, const REAL wDomain,
    const REAL hDomain, TPZVec<int> &nDivs, const MMeshType meshType,
    const bool usingSymmetry)
{
  TPZSimpleTimer timer("create gmesh");
  // dimension of the problem
  constexpr int dim{2};
  // lower left corner of the domain
  TPZManVector<REAL, 3> minX = {0, 0, 0};
  
  // upper right corner of the domain
  TPZManVector<REAL, 3> maxX = {
    usingSymmetry ? wDomain * scale / 2 : wDomain * scale,
    hDomain * scale,
    0};
  
  // whether to create boundary elements
  constexpr bool genBoundEls{true};
  
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,
                                              meshType,genBoundEls);
}

TPZVec<TPZAutoPointer<TPZCompMesh>>
CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
            const TPZVec<int> &matIdVec, const CSTATE ur, const CSTATE er,
            const STATE lambda, const REAL &scale, bool usingSymmetry, ESymType sym)
{
  TPZSimpleTimer timer ("Create cmesh");
  constexpr int dim = 2;
  constexpr bool isComplex{true};


  /*
   First we create the computational mesh associated with the H1 space
   (ez component)*/
  auto * cmeshH1 =new TPZCompMesh(gmesh,isComplex);
  cmeshH1->SetDefaultOrder(pOrder +1);//for deRham compatibility
  cmeshH1->SetDimModel(dim);

  const int volMatId = matIdVec[0];
  //number of state variables in the problem
  constexpr int nState = 1;

  auto dummyMat = new TPZNullMaterial<CSTATE>(volMatId,dim,nState);
  cmeshH1->InsertMaterialObject(dummyMat);

  
  TPZFNMatrix<1, CSTATE> val1(1, 1, 0.);
  TPZManVector<CSTATE,1> val2(1, 0.);
  const int nMats = matIdVec.size();
  TPZBndCond *dummyBC = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == ESymType::PMC ? 1 : 0;
    dummyBC =
      dummyMat->CreateBC(dummyMat, matIdVec[i], bcType, val1, val2);
    cmeshH1->InsertMaterialObject(dummyBC);
  }

  cmeshH1->SetAllCreateFunctionsContinuous();
  cmeshH1->AutoBuild();
  cmeshH1->CleanUpUnconnectedNodes();

  /*
    Then we create the computational mesh associated with the HCurl space
   */
  auto *cmeshHCurl = new TPZCompMesh(gmesh,isComplex);
  cmeshHCurl->SetDefaultOrder(pOrder);
  cmeshHCurl->SetDimModel(dim);
  
  dummyMat = new TPZNullMaterial<CSTATE>(volMatId,dim,nState);
  cmeshHCurl->InsertMaterialObject(dummyMat);

  
  dummyBC = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == ESymType::PMC ? 1 : 0;
    dummyBC =
      dummyMat->CreateBC(dummyMat, matIdVec[i], bcType, val1, val2);
    cmeshHCurl->InsertMaterialObject(dummyBC);
  }

  cmeshHCurl->SetAllCreateFunctionsHCurl();
  cmeshHCurl->AutoBuild();
  cmeshHCurl->CleanUpUnconnectedNodes();

  
  auto *cmeshMF =new TPZCompMesh(gmesh,isComplex);
  
  TPZWaveguideModalAnalysis *matWG  = new TPZWaveguideModalAnalysis(
      volMatId, ur, er, lambda, 1. / scale);
  cmeshMF->InsertMaterialObject(matWG);

  TPZBndCond *bcMat = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == ESymType::PMC ? 1 : 0;
    bcMat =
      matWG->CreateBC(matWG, matIdVec[i], bcType, val1, val2);
    cmeshMF->InsertMaterialObject(bcMat);
  }

  cmeshMF->SetDimModel(dim);
  cmeshMF->SetAllCreateFunctionsMultiphysicElem();

  cmeshMF->AutoBuild();
  cmeshMF->CleanUpUnconnectedNodes();

  TPZManVector<TPZCompMesh*,3> meshVecIn(2);
  meshVecIn[TPZWaveguideModalAnalysis::H1Index()] = cmeshH1;
  meshVecIn[TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl;

  
  TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF);
  TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF);
  TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF);

  cmeshMF->ExpandSolution();
  cmeshMF->ComputeNodElCon();
  cmeshMF->CleanUpUnconnectedNodes();

  TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec(3,nullptr);

  meshVec[0] = cmeshMF;
  meshVec[1 + TPZWaveguideModalAnalysis::H1Index()] = cmeshH1;
  meshVec[1 + TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl;
  return meshVec;
}

void FilterBoundaryEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                             TPZVec<long> &activeEquations, int &neq,
                             int &neqOriginal)
{
  TPZSimpleTimer timer ("Filter dirichlet eqs");
  auto cmesh = meshVec[0];
  TPZManVector<long, 1000> allConnects;
  std::set<long> boundConnects;

  for (int iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if (cel == nullptr) {
      continue;
    }
    if (cel->Reference() == nullptr) {//there is no associated geometric el
      continue;
    }
    TPZBndCond *mat = dynamic_cast<TPZBndCond *>(
        meshVec[0]->MaterialVec()[cel->Reference()->MaterialId()]);
    if (mat && mat->Type() == 0) {//check for dirichlet bcs
      std::set<long> boundConnectsEl;
      std::set<long> depBoundConnectsEl;
      std::set<long> indepBoundConnectsEl;
      cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
      cel->BuildConnectList(boundConnectsEl);
      for (std::set<long>::iterator iT = boundConnectsEl.begin();
           iT != boundConnectsEl.end(); iT++) {
        const long val = *iT;
        if (boundConnects.find(val) == boundConnects.end()) {
          boundConnects.insert(val);
        }
      }
    }
  }
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (boundConnects.find(iCon) == boundConnects.end()) {
      TPZConnect &con = cmesh->ConnectVec()[iCon];
      if (con.HasDependency())
        continue;
      int seqnum = con.SequenceNumber();
      int pos = cmesh->Block().Position(seqnum);
      int blocksize = cmesh->Block().Size(seqnum);
      if (blocksize == 0)
        continue;

      int vs = activeEquations.size();
      activeEquations.Resize(vs + blocksize);
      for (int ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
    }
  }

  neqOriginal = cmesh->NEquations();
  neq = 0;
  auto cmeshHCurl = meshVec[1];
  auto cmeshH1 = meshVec[2];
  int nHCurlEquations = 0, nH1Equations = 0;
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    bool isH1;
    if (boundConnects.find(iCon) == boundConnects.end()) {
      if (cmesh->ConnectVec()[iCon].HasDependency())
        continue;
      int seqnum = cmesh->ConnectVec()[iCon].SequenceNumber();
      int blocksize = cmesh->Block().Size(seqnum);
      if (TPZWaveguideModalAnalysis::H1Index() == 0 && iCon < cmeshH1->NConnects()) {
        isH1 = true;
      } else if (TPZWaveguideModalAnalysis::H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
        isH1 = true;
      } else {
        isH1 = false;
      }
      for (int ieq = 0; ieq < blocksize; ieq++) {
        neq++;
        isH1 == true ? nH1Equations++ : nHCurlEquations++;
      }
    }
  }
  std::cout << "------\tactive eqs\t-------" << std::endl;
  std::cout << "# H1 equations: " << nH1Equations << std::endl;
  std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
  std::cout << "# equations: " << neq << std::endl;
  std::cout << "------\t----------\t-------" << std::endl;
  return;
}