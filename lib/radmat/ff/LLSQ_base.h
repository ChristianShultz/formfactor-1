#ifndef LLSQ_BASE_H_H_GUARD
#define LLSQ_BASE_H_H_GUARD


#include "common_ensemble.h"
#include "semble/semble_vector.h"
#include "adat/handle.h"
#include "formfactor_factory.h"
#include <string>
#include <algorithm>


namespace radmat
{
  // deriving from this polymorphicall will inevitably 
  // give us some freedom later, return handles from the solver routines
  template<typename T>
  struct LLSQRetTypeBase_t
  {
    typedef typename SEMBLE::SembleVector<T> FFType; 
        
    LLSQRetTypeBase_t(const FFType &FFSoln)
      : m_FF(FFSoln)
    {  }

    LLSQRetTypeBase_t(void) 
    : m_FF(SEMBLE::SembleVector<T>())
    {  }

    virtual ~LLSQRetTypeBase_t(void) {}

  public: 
    FFType m_FF;
  };


  // you always get an ensemble of 4 vectors and an ensemble of mat elems
  template<typename T>
  struct LLSQInputType_t
  {
    typedef typename ffGenLLSQSys_t<T>::SemblePInvList_t SemblePInvList_t;
    typedef typename SEMBLE::SembleVector<T> LatticeMatrixElements;
    typedef typename SEMBLE::SembleMatrix<T> KinematicFactorMatrix;

    LLSQInputType_t(const SemblePInvList_t &Momenta,
		    const LatticeMatrixElements  &MatElems,
		    const std::string &matElemID,
		    double precision = 1e-6)
      : m_Momenta(Momenta) , m_MatElems(MatElems) 
    { 
      // sanity -- single precision
      checkQ2(precision);

      // create the linear system 
      ffKinematicFactors_t<T> generator(FormFactorDecompositionFactoryEnv::callFactory(matElemID) );
      m_KFacs = generator(m_Momenta);
    }

    virtual ~LLSQInputType_t(void) {}

  private:
    // sanity
    void checkQ2(double precision)
    {
      typename SemblePInvList_t::const_iterator it;
      itpp::Vec<double> q = mean( m_Momenta.begin()->first - m_Momenta.begin()->second);
      POW2_ASSERT(q.size() == 4);
      const double Q2 =  q[0]*q[0] - (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]) ;
      double mQ2;

      for(it = m_Momenta.begin(); it != m_Momenta.end(); it++)
	{
	  q = mean( it->first - it->second);
	  mQ2 =  q[0]*q[0] - (q[1]*q[1] + q[2]*q[2] + q[3]*q[3]) ;

	  // this is unit normalized precision -- basically just a sanity check
	  if(precision != -1)
	    POW2_ASSERT(fabs(Q2 - mQ2)/fabs(std::max( fabs(Q2), fabs(mQ2) )) < precision);
	}
    }

    // hide ctor
  private:
    LLSQInputType_t(void);

  public:
    SemblePInvList_t m_Momenta;   
    LatticeMatrixElements m_MatElems;
    KinematicFactorMatrix m_KFacs;
  };
  

  template<typename T>
  struct LLSQBaseSolver_t
  {
    typedef typename ADAT::Handle<LLSQRetTypeBase_t<T> > LLSQRetTypeBase_h;
    typedef typename ADAT::Handle<LLSQInputType_t<T> > LLSQInputType_h;
    virtual LLSQRetTypeBase_h operator()(const LLSQInputType_h &) const = 0;
    virtual ~LLSQBaseSolver_t(void) {}
  };


} // close radmat

#endif
