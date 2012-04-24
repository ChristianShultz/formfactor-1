#ifndef LINEAR_SYSTEM_H_H_GUARD
#define LINEAR_SYSTEM_H_H_GUARD

#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_linear_algebra.h"

/**
   @file linear_system.h
   @brief a container to solve Ax = b
 */

namespace linsys
{
  
  /**
     @brief a containet to store the linear system solution
   */
  template<typename T>
  struct LinearSystemSolution
  {
    SembleVector<T> x;
    SembleMatrix<T> cov;
    double residual;
    int nReset;
  }


  /**
     @brief a container to solve Ax = b
   */

  template<typename T>
  struct LinearSystem
  {
    LinearSystem(void); //! not implemented

    //! set up the linear system
    LinearSystem(const SembleVector<T> &b_, const SembleMatrix<T> &A_, const double tolerance_)
      : b(b_) , A(A_) , tolerance(tolerance_)
    {  }

    //! solve the linear system
    LinearSystemSolution<T> solve(void) const
    {
      LinearSystemSolution<T> ret;
      SEMBLE::solveLinearSVD(A,ret.x,b,ret.cov,tolerance,ret.residual,ret.nReset);
      return ret;
    }

  private: // data store
    SembleVector<T> b;
    SembleMatrix<T> A;
    double tolerance;
  }


}

#endif
