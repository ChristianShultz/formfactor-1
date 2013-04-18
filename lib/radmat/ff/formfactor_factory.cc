// formfactor_factory.cc -
//
// Saturday, June  2 2012
//

#include"formfactor_factory.h"
#include <string>
#include <complex>
#include <exception>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "adat/handle.h"

// ffs
#include "lorentzff_PiPi.h"
#include "lorentzff_PiPiStar.h"
#include "lorentzff_PiRho.h"

#include <omp.h>

namespace FacEnv = radmat::FormFactorDecompositionFactoryEnv;
typedef radmat::TheFormFactorDecompositionFactory Factory;


namespace radmat
{

  namespace FormFactorDecompositionFactoryEnv
  {

    namespace 
    {

      // helper function
      template<class T, class U> 
        T* upCast(void)
        {
          T *t = new U();
          POW2_ASSERT(t);
          return t;
        }

      volatile bool registered = false;
    }

    // never played with this toy before so we are just going to 
    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    ADAT::Handle<ffBase_t<std::complex<double> > > callFactory(const std::string &matElemID)
    {
      POW2_ASSERT(FacEnv::registerAll());
      ffBase_t<std::complex<double> > *foo;
      foo = NULL;
      try
      {
        foo = TheFormFactorDecompositionFactory::Instance().createObject(matElemID);
      }
      catch(...)
      {
        POW2_ASSERT(false);
      }

      POW2_ASSERT(foo);
      return ADAT::Handle<ffBase_t<std::complex<double> > >(foo);
    }

    // register the factory "inventory"
    bool registerAll(void)
    {

      bool success = true;
#pragma omp critical
      {

        if(!!!registered)
        {
          success &= Factory::Instance().registerObject(std::string("PiPi"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);
          success &= Factory::Instance().registerObject(std::string("PiPi_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);
          success &= Factory::Instance().registerObject(std::string("PiPiStar_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPiStar::PiPiStar>);
          success &= Factory::Instance().registerObject(std::string("PiRho_0_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<-1> >);
          success &= Factory::Instance().registerObject(std::string("PiRho_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<0> >);
          success &= Factory::Instance().registerObject(std::string("PiRho_0_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<1> >);


          registered = true;
        }
      }
      return success;

    }


  } // close FormFactorDecompositionFactoryEnv

} // close radmat
