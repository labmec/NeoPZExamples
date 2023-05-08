/**
 * @file TPZMatHCurl3D.h
 * @brief Header file for class TPZMatHCurl3D.\n
 */

#ifndef TPZMATHCURL3D_H
#define TPZMATHCURL3D_H


#include "TPZMatBase.h"
#include "TPZMatSingleSpace.h"
/**
 * @ingroup material
 * @brief This class implements the weak statement for a simple problem
 curl (a curl u) - c u = f for illustrating usage of HCurl elements
*/
template<class TVar>
class  TPZMatHCurl3D  :
  public TPZMatBase<TVar,
                    TPZMatSingleSpaceT<TVar>>
{
  using TBase = TPZMatBase<TVar,
                           TPZMatSingleSpaceT<TVar>>;
public:
  /**
     @brief Constructor taking a few material parameters
     @param[in] id Material identifier.
     @param[in] a Coefficient of the curl term
     @param[in] c Coefficient of the mass term
  */
  TPZMatHCurl3D(int id, const TVar a, const TVar c);
  
  TPZMatHCurl3D * NewMaterial() const override;
    
  std::string Name() const override { return "TPZMatHCurl3D"; }

  /** @brief Returns the integrable dimension of the material */
  int Dimension() const override {return 3;}

  [[nodiscard]] int NStateVariables() const override{return 1;}


  /**
     @name SolutionMethods
     @{*/
  //! Variable index of a given solution
  int VariableIndex(const std::string &name) const override;
  //! Number of variables associated with a given solution
  int NSolutionVariables(int var) const override;
  //! Computes the solution at an integration point
  void Solution(const TPZMaterialDataT<TVar> &data,
                int var, TPZVec<TVar> &solout) override;
  /**@}*/
  
  /**
     @name ParamMethods
     @{
  */
  //! Sets the coefficient of the mass termwavelength being analysed
  void SetMassCoeff(const TVar c) {fC = c;}
  //! Gets the current wavelength
  [[nodiscard]] inline TVar GetMassCoeff() const{ return fC;}
  //! Sets the coefficient of the mass termwavelength being analysed
  void SetCurlCoeff(const TVar a) {fA = a;}
  //! Gets the current wavelength
  [[nodiscard]] inline TVar GetCurlCoeff() const{ return fA;}
  /**@}*/
  
  /**
     @name ContributeMethods
     @{
  */
  void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                  TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

  void Contribute(const TPZMaterialDataT<TVar> &data, REAL weight,
                  TPZFMatrix<TVar> &ef) override {}//nothing to be done
  
  void ContributeBC(const TPZMaterialDataT<TVar> &data, REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef,
                    TPZBndCondT<TVar> &bc) override;
  /**@}*/

  [[nodiscard]] int ClassId() const override;
protected:
  TPZMatHCurl3D() = default;
  //! Coefficient of the curl term
  TVar fA{1.};
  //! Coefficient of the mass term
  TVar fC{1.};
};

#endif
