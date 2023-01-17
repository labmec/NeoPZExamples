#include "TPZMatHCurl3D.h"
#include "TPZMaterialDataT.h"
#include <pzaxestools.h>

TPZMatHCurl3D::TPZMatHCurl3D(int id, const STATE c) : TBase(id), fC(c) {}

TPZMatHCurl3D * TPZMatHCurl3D::NewMaterial() const
{
  return new TPZMatHCurl3D(*this);
}



void TPZMatHCurl3D::Contribute(const TPZMaterialDataT<STATE> &data,
                               REAL weight,
                               TPZFMatrix<STATE> &ek,
                               TPZFMatrix<STATE> &ef)
{
  
  const int nshape = data.phi.Rows();
  const auto &phi = data.phi;
  const auto &curl_phi = data.curlphi;
  //evaluates rhs
  constexpr int dim{3};
  TPZManVector<STATE,3> force(dim,0.); 
  if(this->HasForcingFunction()){
    this->fForcingFunction(data.x,force);
  }
  TPZFNMatrix<9,STATE> coeff_mat(3,3,0.), rhs_mat(3,1,0.);
  for(int x = 0; x < 3; x++){
    coeff_mat(x,x) = fC;
    rhs_mat(x,0) = force[x];
  }

  //initializing with correct dimensions
  TPZFNMatrix<3000,STATE> phi_t(3,nshape);
  TPZFNMatrix<3000,STATE> curl_phi_t(nshape,3);

  phi.Transpose(&phi_t);
  curl_phi.Transpose(&curl_phi_t);
    
  ek += ((curl_phi_t * curl_phi) -  phi*(coeff_mat*phi_t))*weight;

  ef += phi*rhs_mat*weight;
}

void TPZMatHCurl3D::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                         REAL weight,
                                         TPZFMatrix<STATE> &ek,
                                         TPZFMatrix<STATE> &ef,
                                         TPZBndCondT<STATE> &bc)
{
  const auto &phi = data.phi;
  const auto& BIG = TPZMaterial::fBigNumber;
    
  const STATE v1 = bc.Val1()(0,0);
  const STATE v2 = bc.Val2()[0];
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
        STATE load{0};
        for(int x = 0; x < 3; x++){
          load += weight * BIG * v2 * phi(i,x);
        }
        ef(i,0) += load;
        for(int j=0;j<nshape;j++){
          STATE stiff{0};
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

int TPZMatHCurl3D::ClassId() const {
  return Hash("TPZMatHCurl3D") ^
    TBase::ClassId() << 1;


}

//! Variable index of a given solution
int TPZMatHCurl3D::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "u") == 0) return 0;
  if( strcmp(name.c_str(), "curl_u") == 0) return 1;
  return TPZMaterial::VariableIndex(name);
}
//! Number of variables associated with a given solution
int TPZMatHCurl3D::NSolutionVariables(int var) const
{
  switch(var){
  case 0: //field (real part)
  case 1: //field (imag val)
    return this->Dimension();
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}
//! Computes the solution at an integration point
void TPZMatHCurl3D::Solution(const TPZMaterialDataT<STATE> &data,
              int var, TPZVec<STATE> &solout)
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