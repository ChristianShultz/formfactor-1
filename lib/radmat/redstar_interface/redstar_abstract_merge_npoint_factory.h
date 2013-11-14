#ifndef REDSTAR_ABSTRACT_MERGE_NPOINT_FACTORY_H
#define REDSTAR_ABSTRACT_MERGE_NPOINT_FACTORY_H 



#include "adat/objfactory.h"
#include "adat/singleton.h"
#include "adat/handle.h"
#include "io/adat_xmlio.h"
#include "redstar_abstract_merge_npoint.h"
#include <string>

namespace radmat
{

  typedef Util::SingletonHolder<
              Util::ObjectFactory<AbsRedstarMergeNPt, 
                                  std::string, 
                                  void,
                                  AbsRedstarMergeNPt* (*)(void), 
                                  Util::StringFactoryError> >
      TheRedstarAbstractMergeNPtFactory; 

  namespace TheRedstarAbstractMergeNPtFactoryEnv
  {
    bool registerAll(void); 
    ADAT::Handle<AbsRedstarMergeNPt> callFactory(const std::string &id); 
  }

} // radmat


#endif /* REDSTAR_ABSTRACT_MERGE_NPOINT_FACTORY_H */
