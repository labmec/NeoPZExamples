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
 curl curl u - c u = for illustrating usage of HCurl elements
*/
class  TPZMatHCurl3D  :
  public TPZMatBase<STATE,TPZMatSingleSpaceT<STATE>>
{
  using TBase = TPZMatBase<STATE,TPZMatSingleSpaceT<STATE>>;
public:
  /**
     @brief Constructor taking a few material parameters
     @param[in] id Material identifier.
     @param[in] c Coefficient of the mass term
  */
  TPZMatHCurl3D(int id, const STATE c);
  
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
  void Solution(const TPZMaterialDataT<STATE> &data,
                int var, TPZVec<STATE> &solout) override;
  /**@}*/
  
  /**
     @name ParamMethods
     @{
  */
  //! Sets the coefficient of the mass termwavelength being analysed
  void SetCoeff(STATE c) {fC = c;}
  //! Gets the current wavelength
  [[nodiscard]] inline STATE GetCoeff() const{ return fC;}
  /**@}*/
  
  /**
     @name ContributeMethods
     @{
  */
  void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;

  void Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                  TPZFMatrix<STATE> &ef) override {}//nothing to be done
  
  void ContributeBC(const TPZMaterialDataT<STATE> &data, REAL weight,
                    TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                    TPZBndCondT<STATE> &bc) override;
  /**@}*/

  [[nodiscard]] int ClassId() const override;
protected:
  TPZMatHCurl3D() = default;
  //! Coefficient of the mass term
  STATE fC{1.};
};

#endif
