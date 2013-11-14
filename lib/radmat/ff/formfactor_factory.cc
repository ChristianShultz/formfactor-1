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

// ffs
#include "lorentzff_PiPi.h"
#include "lorentzff_PiPiStar.h"
#include "lorentzff_PiRho.h"
#include "lorentzff_RhoPi.h"

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
    rHandle<ffBase_t<std::complex<double> > > callFactory(const std::string &matElemID)
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
      return rHandle<ffBase_t<std::complex<double> > >(foo);
    }

    // register the factory "inventory"
    bool registerAll(void)
    {

      bool success = true;

      // guard with critical block for when I inevitably do something stupid
#pragma omp critical
      {

        if(!!!registered)
        {
          // <Pi | jmu | Pi >
          success &= Factory::Instance().registerObject(std::string("PiPi"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);
          success &= Factory::Instance().registerObject(std::string("PiPi_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);
          success &= Factory::Instance().registerObject(std::string("PiPiStar_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPiStar::PiPiStar>);

          // <Pi | jmu | Rho>
          success &= Factory::Instance().registerObject(std::string("PiRho_0_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<-1> >);
          success &= Factory::Instance().registerObject(std::string("PiRho_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<0> >);
          success &= Factory::Instance().registerObject(std::string("PiRho_0_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<1> >);

          // <Rho | jmu | Pi> 
          success &= Factory::Instance().registerObject(std::string("RhoPi_-1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoPi::RhoPi<-1> >);
          success &= Factory::Instance().registerObject(std::string("RhoPi_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoPi::RhoPi<0> >);
          success &= Factory::Instance().registerObject(std::string("RhoPi_1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoPi::RhoPi<1> >);


          registered = true;
        }
      }
      return success;

    }


  } // close FormFactorDecompositionFactoryEnv

} // close radmat
