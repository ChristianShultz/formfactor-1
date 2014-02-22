/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : formfactor_factory.cc

* Purpose :

* Creation Date : 22-02-2014

* Last Modified : Sat 22 Feb 2014 04:25:41 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "formfactor_factory.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_cubic_formfactors.h"
#include <string>
#include <complex>
#include <exception>

namespace FacEnv = radmat::FormFactorDecompositionFactoryEnv;
typedef radmat::TheFormFactorInjectionRecipeFactory Factory;


namespace radmat
{

  namespace FormFactorDecompositionFactoryEnv
  {

    bool registered = false;

    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<FFAbsBase_t > 
      callFactory(const std::string &matElemID)
    {
      FormFactorRecipe_t *foo;
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
      rHandle<FormFactorRecipe_t > recipe(foo);
      if( recipe->id() == ::radmat::Stringify<HelicityFormFactorRecipe_t>() )
      {
        return rHandle<FFAbsBase_t>( new HelicityFormFactor(recipe) ); 
      }
      else
      {
        throw std::string("unknown recipe"); 
        exit(1); 
      }
    }

    // dump all keys in the factory
    std::vector<std::string> 
      all_keys(void)
    {
      return Factory::Instance().keys(); 
    }


    // register the factory "inventory"
    bool registerAll( void )
    {
      bool success = true;

      if(!!!registered)
      {
        success &= ::radmat::HelicityFormFactorDecompositionFactoryEnv::registerAll(); 
        success &= ::radmat::CubicFormFactorDecompositionFactoryEnv::registerAll(); 
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
