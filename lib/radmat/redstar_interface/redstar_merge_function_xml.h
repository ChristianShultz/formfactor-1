#ifndef REDSTAR_MERGE_FUNCTION_XML_H
#define REDSTAR_MERGE_FUNCTION_XML_H 

#include "redstar_abstract_merge_npoint.h"
#include "io/adat_xmlio.h"
#include <string>

namespace radmat
{
  typedef AbstractNamedObject<AbsRedstarMergeNPt>
    AbstractMergeNamedObject;

  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      AbstractMergeNamedObject &);
}

#endif /* REDSTAR_MERGE_FUNCTION_XML_H */
