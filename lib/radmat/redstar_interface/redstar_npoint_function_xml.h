#ifndef REDSTAR_NPOINT_FUNCTION_XML_H
#define REDSTAR_NPOINT_FUNCTION_XML_H 

#include "abstract_named_object.h"
#include "redstar_abstract_block_base.h"
#include "io/adat_xmlio.h"
#include <string> 

namespace radmat
{

  typedef AbstractNamedObject<AbsRedstarXMLInterface_t>
    AbstractBlockNamedObject;

  struct NPointXML
  {
    int version; 
    int N; 
    ADATXML::Array<AbstractBlockNamedObject> npoint;
    std::string ensemble; 
  };

  void read(ADATXML::XMLReader &xml, const std::string &path, NPointXML &npt);
 
  
} // radmat



#endif /* REDSTAR_NPOINT_FUNCTION_XML_H */
