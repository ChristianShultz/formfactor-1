#ifndef REDSTAR_NPOINT_FUNCTION_H
#define REDSTAR_NPOINT_FUNCTION_H 

#include "redstar_abstract_block_base.h"
#include "io/adat_xmlio.h"
#include "adat/handle.h"

namespace radmat
{

  struct AbstractBlockNamedObject
  {
    std::string object_name; 
    ADAT::Handle<AbsRedstarXMLInterface_t> param; 
  };


  struct NPointXML
  {
    int version; 
    int N; 
    ADATXML::Array<AbstractBlockNamedObject> npoint;
  };

  void read(ADATXML::XMLReader &xml, const std::string &path, NPointXML &npt);
  

} // radmat



#endif /* REDSTAR_NPOINT_FUNCTION_H */
