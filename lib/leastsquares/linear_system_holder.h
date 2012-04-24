#ifndef LINEAR_SYSTEM_HOLDER_H_H_GUARD
#define LINEAR_SYSTEM_HOLDER_H_H_GUARD

#include "kinematic_factors/kinematic_factors.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"

/**
   @file linear_system_holder.h
   @brief hold the matrix of kinematic factors
 */

namespace linsys
{
  
  /**
     @brief a container to objectify some data and make it easier to pass 
   */
  struct LorentzIndices
  {
    LorentzIndicies(void)
    : zero(false) , one(false) , two(false) , three(false)
    {  }

    LorentzIndicies(const bool a, const bool b, const bool c, const bool d)
      : zero(a) , one(b) , two(c) , three(d)
    {  }

  public:
    bool zero;
    bool one; 
    bool two;
    bool three;
  }


  /**
     @brief hold the ensemble matrix of non-zero kinematic factors
   */
  template<typename T>
  struct LinSysHolder
  {
    //! constructor
    LinSysHolder(void)
    : init(false)
    {  }
    
    //! add a set of rows to the linear system
    void add(const LorentzIndicies &l, const EnsembleRowHolder &rows)
    {
      if(!!!init)
	{
	  init = true;
	  linsys_ = SembleMatrix<T>(rows[0].bins(),0,rows[0].getN()); 
	}

      if(l.zero)
	linsys_.append_row(rows[0]);
      if(l.one)
	linsys_.append_row(rows[1]);
      if(l.two)
	linsys_.append_row(rows[2]);
      if(l.three)
	linsys_.append_row(rows[3]);
    }
    
    //! return the linear system
    SEMBLE::SembleMatrix<T> getLinSys(void)
    {
      return linsys_;
    }

  private: // data store
    bool init;
    SEMBLE::SembleMatrix<T> linsys_;
  }


}

#endif
