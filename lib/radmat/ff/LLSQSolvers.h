#ifndef LLSQSOLVERS_H_H_GUARD
#define LLSQSOLVERS_H_H_GUARD

#include "LLSQ_base.h"
#include <string>
#include <complex>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "adat/handle.h"
#include "semble/semble_matrix.h"
#include "semble/semble_linear_algebra.h"

namespace radmat
{

  // factory
  typedef Util::SingletonHolder<
    Util::ObjectFactory<LLSQBaseSolver_t<std::complex<double> >,
			std::string,
			void,
			LLSQBaseSolver_t<std::complex<double> >* (*)(),
			Util::StringFactoryError > >
  TheLLSQSolverFactory;

  namespace LLSQSolverFactoryEnv
  {
    bool registerAll();
    ADAT::Handle<LLSQBaseSolver_t<std::complex<double> > > callFactory(const std::string &solnID);
  }


  /////////////////////////////////////////////////////////////////////////////

  // solver functors
  template<typename T> 
  struct LLSQSolverLU_t : public LLSQBaseSolver_t<T>
  {
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h & input) const
    {
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);
      SEMBLE::SembleMatrix<T> Kinv;
      SEMBLE::inv(input->m_KFacs,Kinv);
      foo->m_FF = Kinv*(input->m_MatElems);
      return LLSQRetTypeBase_h(foo);
    }
  };



}

#endif
