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
      
      bool registered = false;
    }

    // never played with this toy before so we are just going to 
    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    ADAT::Handle<ffBase_t<std::complex<double> > > callFactory(const std::string &matElemID)
    {
      POW2_ASSERT(FacEnv::registerAll());
      ffBase_t<std::complex<double> > *foo;

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

      if(!!!registered)
	{
	  success &= Factory::Instance().registerObject(std::string("PiPi"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);


	  registered = true;
	}

      return success;
    }


  } // close FormFactorDecompositionFactoryEnv

} // close radmat
