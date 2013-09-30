#ifndef BUILD_CORRELATORS_XML_H
#define BUILD_CORRELATORS_XML_H



#include "g_parity_world_xml.h"
#include "g_parity_world.h"
#include "g_parity_world_generate_redstar_xml.h"
#include "radmat_database_interface.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include <string>
#include <iostream>
#include <iomanip>



namespace radmat
{
  struct ThreePointCorrXMLIni_t
  {
    struct RenormalizationProp
    {
      double RGE_prop;  // flavor space wick contraction normalization 
      double Z_t; 
      double Z_s; 
    };

    gParityWorld::GParityContinuumMatElemXML continuumMatElemXML;
    std::string source_id;
    std::string sink_id;
    bool gParitySymmetry;
    bool cubicSymmetry;
    double maSource;
    double maSink;  

    RenormalizationProp renormalization; 
  };

  //! write a renormalization prop to a string 
  std::string toString(const ThreePointCorrXMLIni_t::RenormalizationProp &);

  //! stream it
  std::ostream& operator<<(std::ostream&, const ThreePointCorrXMLIni_t::RenormalizationProp &); 

  //! read it
  void read(ADATXML::XMLReader &xml, 
      const std::string &path, 
      ThreePointCorrXMLIni_t::RenormalizationProp &); 

  //! write it
  void write(ADATXML::XMLWriter &xml, 
      const std::string &path, 
      const ThreePointCorrXMLIni_t::RenormalizationProp &);

  //! write it to a string
  std::string toString(const ThreePointCorrXMLIni_t &);

  //! stream it
  std::ostream& operator<<(std::ostream&, const ThreePointCorrXMLIni_t &);

  //! xml reader
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrXMLIni_t &);

  //! xml writer
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrXMLIni_t &);


  struct ThreePointCorrIni_t
  {
    ThreePointCorrXMLIni_t threePointCorrXMLIni;
    radmatDBProp_t radmatDBProp;
    std::string matElemID;
    double xi;
    int L_s;   
  };

  // boiler plate stuff
  std::string toString(const ThreePointCorrIni_t &);
  std::ostream& operator<<(std::ostream&, const ThreePointCorrIni_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrIni_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrIni_t &);

 
} // radmat 

#endif /* BUILD_CORRELATORS_XML_H */
