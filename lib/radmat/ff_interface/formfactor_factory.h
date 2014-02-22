#ifndef FORMFACTOR_FACTORY_H
#define FORMFACTOR_FACTORY_H 


#include "formfactor.h"
#include <string>
#include <complex>
#include <vector>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/handle.h"

//
//
//    helicities are a direct map back to cont_spin ffs
//    cubics are weighted sums of cont_spin ffs
//
//


namespace radmat
{

  // inject recipes upon instantiation 
  typedef Util::SingletonHolder<
    Util::ObjectFactory<FormFactorRecipe_t,
			std::string,
			void,
			FormFactorRecipe_t* (*)(void),
			Util::StringFactoryError> >
  TheFormFactorInjectionRecipeFactory;


  // return a recipe injected form factor
  namespace FormFactorDecompositionFactoryEnv
  {
    bool registerAll( void );
    rHandle<FFAbsBase_t> callFactory(const std::string &matElemID);
    std::vector<std::string> all_keys(void); 
  }

}



#endif /* FORMFACTOR_FACTORY_H */
