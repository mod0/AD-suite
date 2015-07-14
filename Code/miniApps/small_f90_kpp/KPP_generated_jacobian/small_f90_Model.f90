MODULE small_f90_Model

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Completely defines the model small_f90
!    by using all the associated modules
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE small_f90_Precision
  USE small_f90_Parameters
  USE small_f90_Global
  USE small_f90_Function
  USE small_f90_Integrator
  USE small_f90_Rates
  USE small_f90_Jacobian
  USE small_f90_Hessian
  USE small_f90_Stoichiom
  USE small_f90_LinearAlgebra
  USE small_f90_Monitor
  USE small_f90_Util

END MODULE small_f90_Model

