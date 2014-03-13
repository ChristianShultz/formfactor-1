#ifndef FORMFACTOR_INVARIANTS_H
#define FORMFACTOR_INVARIANTS_H 

#include <string>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"

namespace radmat
{
  struct FFRep_p; 
  REGISTER_STRINGIFY_TYPE( FFRep_p ); 


  struct FFRep_p
  {
    virtual ~FFRep_p() {}
    virtual std::string rep_type(void) const {return Stringify<FFRep_p>();}
    virtual std::string rep_id(void) const = 0; 
    virtual int rep_row(void) const = 0; 
  };


  // factory to hold them all
  typedef Util::SingletonHolder<
    Util::ObjectFactory<FFRep_p,
			std::string,
			void,
			FFRep_p* (*)(void),
			Util::StringFactoryError> >
  TheFormFactorInvariantsFactory;

  namespace FormFactorInvariantsFactoryEnv
  {
    bool registerAll(void);

    // nb keys are class names
    rHandle<FFRep_p> callFactory(const std::string &id);
  }

}; 


#endif /* FORMFACTOR_INVARIANTS_H */
