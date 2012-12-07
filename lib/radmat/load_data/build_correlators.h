#ifndef BUILD_CORRELATORS_H
#define BUILD_CORRELATORS_H


#include "simple_world.h"
#include "radmat_database_interface.h"
#include "generate_redstar_xml.h"  
#include "radmat/llsq/llsq_q2_pack.h"
#include "semble/semble_key_val_db.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include <string>
#include <iostream>


namespace radmat
{


  struct ThreePointCorrXMLIni_t
  {
    simpleWorld::ContinuumMatElemXML continuumMatElemXML;
    std::string source_id;
    std::string sink_id;
    bool isDiagonal; 
  };


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
    double L_s;   
  };

  // boiler plate stuff
  std::string toString(const ThreePointCorrIni_t &);
  std::ostream& operator<<(std::ostream&, const ThreePointCorrIni_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrIni_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrIni_t &);


  struct BuildCorrelators
  {
    BuildCorrelators(void) : have_ini(false) {}
    BuildCorrelators(const ThreePointCorrIni_t &ini) : have_ini(true) , m_ini(ini) {}

    void load(const ThreePointCorrIni_t &ini)
    {
      have_ini = true;
      m_ini = ini; 
    }

    std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > build_correlators(const ThreePointCorrIni_t &ini)
    {
      load(ini);
      return build_correlators(); 
    }
    std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > build_correlators(void); 

    bool have_ini;
    ThreePointCorrIni_t m_ini;
  };

};








#endif /* BUILD_CORRELATORS_H */
