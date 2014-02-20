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
#include "lorentzff_canonical_PiRho.h"
#include "lorentzff_canonical_RhoPi.h"
#include "lorentzff_canonical_RhoRho.h"

#include <omp.h>

namespace FacEnv = radmat::FormFactorDecompositionFactoryEnv;
typedef radmat::TheFormFactorDecompositionFactory Factory;


namespace radmat
{

  namespace FormFactorDecompositionFactoryEnv
  {

    // helper function
    template<class T, class U> 
      T* upCast(void)
      {
        T *t = new U();
        POW2_ASSERT(t);
        return t;
      }

    template<typename T> 
      bool 
      do_reg(const std::string &reg_id, T* (*ptr)())
      {
        bool reg = Factory::Instance().registerObject(reg_id,ptr); 

        if ( !!! reg ) 
        {
          std::cout << __PRETTY_FUNCTION__ << ": reg error for " << reg_id << std::endl;
        } 
        return reg; 
      }

    bool registered = false;

    // never played with this toy before so we are just going to 
    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<FFAbsBase_t > callFactory(const std::string &matElemID)
    {
      FFAbsBase_t *foo;
      foo = NULL;
      try
      {
        foo = TheFormFactorDecompositionFactory::Instance().createObject(matElemID);
      }
      catch(std::exception &e)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << e.what(); 
        throw e; 
      }
      catch(std::string &s)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << s << std::endl;
        throw s;
      }
      catch(...)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": some error" << std::endl;
        POW2_ASSERT(false); 
      }

      POW2_ASSERT(foo);
      return rHandle<FFAbsBase_t >(foo);
    }

    // register the factory "inventory"
    bool registerAll(void)
    {

      bool success = true;

      if(!!!registered)
      {
        // <Pi | jmu | Pi >
        // success &= Factory::Instance().registerObject(std::string("PiPi"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);

        success &= do_reg(std::string("PiPi"), FacEnv::upCast<FFAbsBase_t ,radmat::PiPi::PiPi>);
        success &= do_reg(std::string("PiPiStar"),FacEnv::upCast<FFAbsBase_t ,radmat::PiPiStar::PiPiStar>);

        // <Pi | jmu | Rho>
        success &= do_reg(std::string("CanonicalPiRho"),FacEnv::upCast<FFAbsBase_t, radmat::CanonicalPiRho::CanonicalPiRho >);

        // <Rho | jmu | Pi> 
        success &= do_reg(std::string("CanonicalRhoPi"),FacEnv::upCast<FFAbsBase_t, radmat::CanonicalRhoPi::CanonicalRhoPi>);

        // <Rho | jum | Rho> 
        success &= do_reg(std::string("RhoRho"),FacEnv::upCast<FFAbsBase_t, radmat::RhoRho::RhoRho >);

        registered = true;
      }

      if( !!! success )
      {
        throw std::string("failed to reg in FormFactorDecompositionFactoryEnv"); 
      }

      return success;

    }


  } // close FormFactorDecompositionFactoryEnv

} // close radmat
