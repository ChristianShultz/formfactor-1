// LLSQSolvers.cc -
//
// Saturday, June  2 2012
//

#include"llsq_solvers.h"
#include <string>
#include <complex>
#include <exception>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "adat/handle.h"
#include "radmat/utils/pow2assert.h"

#include <omp.h>

namespace FacEnv = radmat::LLSQSolverFactoryEnv;
typedef radmat::TheLLSQSolverFactory Factory;

namespace radmat
{

  namespace LLSQSolverFactoryEnv
  {

    namespace
    {
      template<class T, class U>
        T* upCast(void)
        {
          T *t = new U();
          POW2_ASSERT(t);
          return t;
        }

      volatile bool registered = false;
    } // close anonymous 



    bool registerAll(void)
    {

      bool success = true;
#pragma omp critical
      {

        if(!!!registered)
        {
          success &= Factory::Instance().registerObject(std::string("LU"),
              FacEnv::upCast<LLSQBaseSolver_t<std::complex<double> >, 
              LLSQSolverLU_t<std::complex<double> > > );
          success &= Factory::Instance().registerObject(std::string("SVDMakeSquare"),
              FacEnv::upCast<LLSQBaseSolver_t<std::complex<double> >, 
              LLSQSolverSVDMakeSquare_t<std::complex<double> > > );

          registered = true;
        }
      } // critical

      return success;
    }


    // interface between the factory and the outside world for ease of use
    ADAT::Handle<LLSQBaseSolver_t<std::complex<double> > > callFactory(const std::string &solnID)
    {
      POW2_ASSERT(FacEnv::registerAll());
      LLSQBaseSolver_t<std::complex<double> > *foo;
      foo = NULL;
      try
      {
        foo = Factory::Instance().createObject(solnID);
      }
      catch(...)
      {
        POW2_ASSERT(false);
      }

      POW2_ASSERT(foo);
      return ADAT::Handle<LLSQBaseSolver_t<std::complex<double> > >(foo);
    }


  } // Env

} // radmat
