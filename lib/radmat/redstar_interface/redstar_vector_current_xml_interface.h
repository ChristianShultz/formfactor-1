#ifndef REDSTAR_VECTOR_CURRENT_XML_INTERFACE_H
#define REDSTAR_VECTOR_CURRENT_XML_INTERFACE_H 

#include "redstar_abstract_xml_interface.h"


  //
  //
  // A BASIC gamma^mu VECTOR CURRENT XML INTERFACE
  //
  //

namespace radmat
{



  struct RedstarVectorCurrentXML;
  REGISTER_STRINGIFY_TYPE(RedstarVectorCurrentXML);


  struct RedstarVectorCurrentXML 
    : public AbsRedstarXMLInterface_t 
  {

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarVectorCurrentXML>());
    }

    virtual std::string write(void) const; 

    // a photon fragment 
    struct pfrag
    {
      double coeff_r; 
      double coeff_i; 
      std::string name; 
    };

    struct insertion
    {
      bool active; 
      bool smearedP; 
      ADATXML::Array<pfrag> photons;  
    };

    int pmin; 
    int pmax; 
    int t_slice; 
    insertion time; 
    insertion space; 
  };

}

#endif /* REDSTAR_VECTOR_CURRENT_XML_INTERFACE_H */
