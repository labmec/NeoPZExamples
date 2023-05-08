#include "TPZMatHCurl3D.h"
#include "TPZMaterialDataT.h"
#include <pzaxestools.h>


template<class TVar>
TPZMatHCurl3D<TVar>::TPZMatHCurl3D(int id, const TVar a, const TVar c) :
  TBase(id), fA(a), fC(c) {}

template<class TVar>
TPZMatHCurl3D<TVar> * TPZMatHCurl3D<TVar>::NewMaterial() const
{
  return new TPZMatHCurl3D(*this);
}


template<class TVar>
void TPZMatHCurl3D<TVar>::Contribute(const TPZMaterialDataT<TVar> &data,
                               REAL weight,
                               TPZFMatrix<TVar> &ek,
                               TPZFMatrix<TVar> &ef)
{
  
  const int nshape = data.phi.Rows();
  
  //evaluates rhs
  constexpr int dim{3};
  TPZManVector<TVar,3> force(dim,0.); 
  if(this->HasForcingFunction()){
    this->fForcingFunction(data.x,force);
  }
  TPZFNMatrix<9,TVar> a_mat(3,3,0.), c_mat(3,3,0.), rhs_mat(3,1,0.);
  for(int x = 0; x < 3; x++){
    a_mat(x,x) = fA;
    c_mat(x,x) = fC;
    rhs_mat(x,0) = force[x];
  }

  TPZFNMatrix<200,TVar> phi(nshape,3,0.), curl_phi(3,nshape,0.);
  for(int i = 0; i < nshape; i++){
    for(int x = 0; x < 3; x++){
      phi(i,x) = data.phi.GetVal(i,x);
      curl_phi(x,i) = data.curlphi.GetVal(x,i);
    }
  }

  
  //initializing with correct dimensions
  TPZFNMatrix<3000,TVar> phi_t(3,nshape);
  TPZFNMatrix<3000,TVar> curl_phi_t(nshape,3);

  phi.Transpose(&phi_t);
  curl_phi.Transpose(&curl_phi_t);
    
  ek += ((curl_phi_t * (a_mat * curl_phi)) -  phi*(c_mat*phi_t))*weight;

  ef += phi*rhs_mat*weight;
}

template<class TVar>
void TPZMatHCurl3D<TVar>::ContributeBC(const TPZMaterialDataT<TVar> &data,
                                       REAL weight,
                                       TPZFMatrix<TVar> &ek,
                                       TPZFMatrix<TVar> &ef,
                                       TPZBndCondT<TVar> &bc)
{
  const auto &phi = data.phi;
  const auto& BIG = TPZMaterial::fBigNumber;
    
  const TVar v1 = bc.Val1()(0,0);
  const TVar v2 = bc.Val2()[0];
  constexpr STATE tol = std::numeric_limits<STATE>::epsilon();
  if(std::abs(v2) > tol){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nThis method supports only homogeneous boundary conditions.\n";
    std::cout<<"Stopping now..."<<std::endl;
    DebugStop();
  }
  switch ( bc.Type() )
    {
    case 0:{
      const int nshape=phi.Rows();
      for(int i = 0 ; i<nshape ; i++){
        TVar load{0};
        for(int x = 0; x < 3; x++){
          load += weight * BIG * v2 * phi(i,x);
        }
        ef(i,0) += load;
        for(int j=0;j<nshape;j++){
          TVar stiff{0};
          for(int x = 0; x < 3; x++){
            stiff += phi(i,x) * phi(j,x) * BIG ;
          }
          ek(i,j) += stiff*weight;
        }
      }
      break;
    }
    case 1:
      ///PMC condition just adds zero to both matrices. nothing to do here....
      break;
    default:
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
      PZError<<"Stopping now..."<<std::endl;
      DebugStop();
      break;
    }
}

template<class TVar>
int TPZMatHCurl3D<TVar>::ClassId() const {
  return Hash("TPZMatHCurl3D") ^
    TBase::ClassId() << 1;


}

template<class TVar>
int TPZMatHCurl3D<TVar>::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "u") == 0) return 0;
  if( strcmp(name.c_str(), "curl_u") == 0) return 1;
  return TPZMaterial::VariableIndex(name);
}

template<class TVar>
int TPZMatHCurl3D<TVar>::NSolutionVariables(int var) const
{
  switch(var){
  case 0: //field (real part)
  case 1: //field (imag val)
    return this->Dimension();
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}

template<class TVar>
void TPZMatHCurl3D<TVar>::Solution(const TPZMaterialDataT<TVar> &data,
              int var, TPZVec<TVar> &solout)
{

  const auto &sol = data.sol[0];
  const auto &curlsol = data.curlsol[0];
  switch (var){
  case 0:
    for(auto x = 0; x < 3; x++){solout[x] = sol[x];}
    break;
  case 1:
    for(auto x = 0; x < 3; x++){solout[x] = curlsol[x];}
    break;
  }
}

template class TPZMatHCurl3D<STATE>;
template class TPZMatHCurl3D<CSTATE>;