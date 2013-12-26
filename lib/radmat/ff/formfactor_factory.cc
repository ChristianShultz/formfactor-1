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
    rHandle<ffBase_t<std::complex<double> > > callFactory(const std::string &matElemID)
    {
      ffBase_t<std::complex<double> > *foo;
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
      return rHandle<ffBase_t<std::complex<double> > >(foo);
    }

    // register the factory "inventory"
    bool registerAll(void)
    {

      bool success = true;

      if(!!!registered)
      {
        // <Pi | jmu | Pi >
        // success &= Factory::Instance().registerObject(std::string("PiPi"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);

        success &= do_reg(std::string("PiPi"), FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);
        success &= do_reg(std::string("PiPi_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPi::PiPi>);
        success &= do_reg(std::string("PiPiStar_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> > ,radmat::PiPiStar::PiPiStar>);

        // <Pi | jmu | Rho>
        success &= do_reg(std::string("PiRho_0_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<-1> >);
        success &= do_reg(std::string("PiRho_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<0> >);
        success &= do_reg(std::string("PiRho_0_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::PiRho::PiRho<1> >);
        // <Pi | jmu | Rho>
        success &= do_reg(std::string("CanonicalPiRho_0_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::CanonicalPiRho::CanonicalPiRho<-1> >);
        success &= do_reg(std::string("CanonicalPiRho_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::CanonicalPiRho::CanonicalPiRho<0> >);
        success &= do_reg(std::string("CanonicalPiRho_0_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::CanonicalPiRho::CanonicalPiRho<1> >);


        // <Rho | jmu | Pi> 
        success &= do_reg(std::string("RhoPi_-1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoPi::RhoPi<-1> >);
        success &= do_reg(std::string("RhoPi_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoPi::RhoPi<0> >);
        success &= do_reg(std::string("RhoPi_1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoPi::RhoPi<1> >);
        // <Rho | jmu | Pi> 
        success &= do_reg(std::string("CanonicalRhoPi_-1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::CanonicalRhoPi::CanonicalRhoPi<-1> >);
        success &= do_reg(std::string("CanonicalRhoPi_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::CanonicalRhoPi::CanonicalRhoPi<0> >);
        success &= do_reg(std::string("CanonicalRhoPi_1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::CanonicalRhoPi::CanonicalRhoPi<1> >);

        // <Rho | jum | Rho> 
        success &= do_reg(std::string("RhoRho_1_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<1,-1> >);
        success &= do_reg(std::string("RhoRho_1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<1,0> >);
        success &= do_reg(std::string("RhoRho_1_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<1,1> >);

        success &= do_reg(std::string("RhoRho_0_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<0,-1> >);
        success &= do_reg(std::string("RhoRho_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<0,0> >);
        success &= do_reg(std::string("RhoRho_0_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<0,1> >);

        success &= do_reg(std::string("RhoRho_-1_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<1,1> >);
        success &= do_reg(std::string("RhoRho_-1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<-1,0> >);
        success &= do_reg(std::string("RhoRho_-1_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::RhoRho::RhoRho<-1,1> >);

        // <Rho | jum | Rho> 
        success &= do_reg(std::string("VecVec_1_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<1,-1> >);
        success &= do_reg(std::string("VecVec_1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<1,0> >);
        success &= do_reg(std::string("VecVec_1_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<1,1> >);

        success &= do_reg(std::string("VecVec_0_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<0,-1> >);
        success &= do_reg(std::string("VecVec_0_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<0,0> >);
        success &= do_reg(std::string("VecVec_0_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<0,1> >);

        success &= do_reg(std::string("VecVec_-1_-1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<1,1> >);
        success &= do_reg(std::string("VecVec_-1_0"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<-1,0> >);
        success &= do_reg(std::string("VecVec_-1_1"),FacEnv::upCast<ffBase_t<std::complex<double> >, radmat::VecVec::VecVec<-1,1> >);

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
