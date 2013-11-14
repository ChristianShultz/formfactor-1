#ifndef FORMFACTOR_FACTORY_H_H_GUARD
#define FORMFACTOR_FACTORY_H_H_GUARD

#include "ff_gen_llsq_row.h"
#include <string>
#include <complex>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/handle.h"

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
    rHandle<ffBase_t<std::complex<double> > > callFactory(const std::string &matElemID);
  }

}

#endif
