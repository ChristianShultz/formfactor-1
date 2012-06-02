#ifndef FORMFACTOR_FACTORY_H_H_GUARD
#define FORMFACTOR_FACTORY_H_H_GUARD

#include "common_cfg.h"
#include "common_ensemble.h"
#include <string>
#include <complex>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "adat/handle.h"

namespace radmat
{

  typedef Util::SingletonHolder<
    Util::ObjectFactory<ffBase_t<std::complex<double> >,
			std::string,
			void,
			ffBase_t<std::complex<double> >* (*)(void),
			Util::StringFactoryError> >
  TheFormFactorDecompositionFactory;


  namespace FormFactorDecompositionFactoryEnv
  {
    bool registerAll();
    ADAT::Handle<ffBase_t<std::complex<double> > > callFactory(const std::string &matElemID);
  }

}

#endif
