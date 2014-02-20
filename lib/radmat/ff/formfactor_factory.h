#ifndef FORMFACTOR_FACTORY_H_H_GUARD
#define FORMFACTOR_FACTORY_H_H_GUARD

#include "formfactor_abs_base_cfg.h"
#include <string>
#include <complex>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/handle.h"

namespace radmat
{

  enum FFMODE
  {
    HELICITY, 
    CUBIC
  };

  typedef Util::SingletonHolder<
    Util::ObjectFactory<ffBase_t<std::complex<double> >,
			std::string,
			void,
			FFAbsBase_t* (*)(void),
			Util::StringFactoryError> >
  TheFormFactorDecompositionFactory;


  namespace FormFactorDecompositionFactoryEnv
  {
    bool registerAll( const FFMODE );
    rHandle<FFAbsBase_t> callFactory(const std::string &matElemID);
  }

}

#endif
