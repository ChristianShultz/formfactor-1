#ifndef COMMON_ENSEMBLE_H_H_GUARD
#define COMMON_ENSEMBLE_H_H_GUARD

#include "common_cfg.h"
#include "adat/handle.h"
#include "semble/semble_matrix.h"
#include "semble/semble_vector.h"
#include "radmat/utils/pow2assert.h"

/**
   @file common_ensemble.h
   @brief set up the linear least squares system for inversion
 */

namespace radmat
{

  // save some typing
  typedef std::pair<SEMBLE::SembleVector<double> , 
		    SEMBLE::SembleVector<double> > SemblePInv_t;


  /**
     @brief generate the kinematic factor matrix for one measurement,
     
     @details ie: a 4 x n form factors matrix where the row index is the lorentz index.
          These then need to get stitched together higher up 
     
   */
  template<typename T>
  struct ffKinematicFactors_t
  {
    // save some typing
    typedef typename ADAT::Handle<ffBase_t<T> > ffBase_h;
    typedef typename SEMBLE::SembleMatrix<T> KinematicFactorMatrix;

    // the only available constructor
    ffKinematicFactors_t(ffBase_h &KFacGen)
    : m_KFacGen(KFacGen) 
    {  }

    ~ffKinematicFactors_t(void) {} // handle cleans itself up

    // basically generate the 4 X (n multipole) matrix of kinematic factors, 
    // the row index is the lorentz index of the lattice matrix element
    KinematicFactorMatrix genFactors(const SemblePInv_t &moms)
    {
      SEMBLE::SembleVector<double> p_f(moms.first), p_i(moms.second);
      int nfacs = m_KFacGen->nFacs();
      int nbins = p_f.getB();

      POW2_ASSERT_DEBUG( (nbins == p_i.getB()) && (nfacs > 0) );

      KinematicFactorMatrix KF(nbins,4,nfacs);
      KF.zeros();

      // scale down
      p_f.rescaleSembleDown();
      p_i.rescaleSembleDown();

      for(int bin = 0; bin < nbins; bin++)
	KF[bin] = (*m_KFacGen)(toTensor<double>(p_f[bin]),toTensor<double>(p_i[bin]));

      ////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////
      /*
	TO DO -- put the isZero()? type check in here before the return statement
      */
      ////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////

      // scale up
      KF.rescaleSembleUp();
    
      return KF;
    }

    // hide ctor    
  private:
    ffKinematicFactors_t(void);
    ffKinematicFactors_t(const ffKinematicFactors_t &o);
    ffKinematicFactors_t& operator=(const ffKinematicFactors_t<T> &o);

  private:
    ffBase_h m_KFacGen;
  };

  /**
     @brief stitch together all of the individual lorentz pieces of the linear system
     @details the output of this is the actual matrix used in inverting K to solve K^-1 M = FF
   */

  template<typename T>
  struct ffGenLLSQSys_t
  {
    // save some typing
    typedef typename ffKinematicFactors_t<T>::ffBase_h ffBase_h;
    typedef typename ffKinematicFactors_t<T>::KinematicFactorMatrix KinematicFactorPiece; 
    typedef std::list<SemblePInv_t> SemblePInvList_t;

    // the only available constructor -- KinematicFactorGenerator is polymorphic and specifies the 
    // type of matrix element we are decomposing ie:  PiPi  <-----> <0^+ | j_mu | 0^+>  
    ffGenLLSQSys_t(const ffBase_h &KinematicFactorGenerator)
    : m_KFacGen(KinematicFactorGenerator)
    {  }

    // this is really something of a functor class -- this constructs the matrix K in the linear system
    typename ffKinematicFactors_t<T>::KinematicFactorMatrix operator()(const SemblePInvList_t &list)
    {
      SemblePInvList_t::const_iterator it = list.begin();

      KinematicFactorPiece st;
      st = (*m_KFacGen)(*it);
      
      it++;

      while(it != list.end())
	{
	  st.append_row( (*m_KFacGen)(*it) );
	  it++;
	}

      return st;
    }
    
    // hide ctor
  private:
    ffGenLLSQSys_t(void);
    ffGenLLSQSys_t(const ffGenLLSQSys_t &o);
    ffGenLLSQSys_t& operator=(const ffGenLLSQSys_t<T> &o);
    
  private:
    ffKinematicFactors_t<T> m_KFacGen;
  };



}

#endif
