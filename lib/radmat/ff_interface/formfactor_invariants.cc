/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : formfactor_invariants.cc

* Purpose :

* Creation Date : 11-03-2014

* Last Modified : Wed 12 Mar 2014 11:28:18 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "formfactor_invariants.h"
#include "formfactor_spherical_invariants.h"
#include "formfactor_cubic_invariants.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"



namespace FacEnv = radmat::FormFactorInvariantsFactoryEnv;
typedef radmat::TheFormFactorInvariantsFactory Factory;


namespace radmat
{

  namespace FormFactorInvariantsFactoryEnv
  {
    namespace 
    {
      bool local_registration = false; 
    }

    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 
      if( !!! local_registration )
      {
        success &= SpherInvariantsFactoryEnv::registerAll(); 
        success &= CubicInvariantsFactoryEnv::registerAll(); 
        local_registration = true; 
      }

      return success; 
    }

    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<FFRep_p > callFactory(const std::string &matElemID)
    {
      FFRep_p *foo;
      foo = NULL;
      try
      {
        foo = Factory::Instance().createObject(matElemID);
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

      // not a null pointer
      POW2_ASSERT(foo);
      return rHandle<FFRep_p >(foo);
    }


  }// Fac Env

}
