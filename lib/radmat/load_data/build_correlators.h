#ifndef BUILD_CORRELATORS_H
#define BUILD_CORRELATORS_H


#include "simple_world.h"
#include "radmat_database_interface.h"
#include "generate_redstar_xml.h"  
#include "radmat/llsq/llsq_q2_pack.h"
#include "radmat/llsq/llsq_multi_data.h"
#include "radmat_overlap_key_val_db.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include <string>
#include <iostream>
#include <iomanip>

namespace radmat
{


  struct ThreePointCorrXMLIni_t
  {
    simpleWorld::ContinuumMatElemXML continuumMatElemXML;
    std::string source_id;
    std::string sink_id;
    bool isDiagonal;
    bool isProjected;
    double maSource;
    double maSink;  
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
    int L_s;   
  };

  // boiler plate stuff
  std::string toString(const ThreePointCorrIni_t &);
  std::ostream& operator<<(std::ostream&, const ThreePointCorrIni_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrIni_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrIni_t &);



  struct LatticeMultiDataTag
  {

    LatticeMultiDataTag(void)
    {
      E_f.resize(1); 
      E_f = SEMBLE::toScalar(double(0.));
      E_i = E_f; 
    }


    ENSEM::EnsemReal Q2(void) const
    {
      double pp(0); 
      pp = mom_fac*mom_fac*((p_f[0] - p_i[0])*(p_f[0] - p_i[0])
          + (p_f[1] - p_i[1])*(p_f[1] - p_i[1])
          + (p_f[2] - p_i[2])*(p_f[2] - p_i[2]));


      return ( - (E_f-E_i)*(E_f-E_i) + SEMBLE::toScalar(pp));
    }

    void print_me(void) const
    {
      std::cout << file_id << " " << jmu << " " << mat_elem_id << std::endl;  
    }


    std::string mom_string(void) const
    {
      std::stringstream ss;
      ss << "pf = " << p_f[0] << "," << p_f[1] << ","
        << p_f[2] << "  pi = "  << p_i[0] << "," 
        << p_i[1] << "," << p_i[2] ;
      return ss.str();
    }


    std::string E_string(void) const
    {
      std::stringstream ss;
      ss << "E_f = " << std::setw(3) << SEMBLE::toScalar(ENSEM::mean(E_f)) 
        << " E_i = " << std::setw(3) << SEMBLE::toScalar(ENSEM::mean(E_i));
      return ss.str(); 
    }


    // tags

    std::string file_id; // some unique string telling us what this is

    // for llsq system
    int jmu;
    std::string mat_elem_id; 
    Array<int> p_f;
    Array<int> p_i;
    ENSEM::EnsemReal E_f;
    ENSEM::EnsemReal E_i;
    double mom_fac; 
  };


  typedef LLSQMultiData<LatticeMultiDataTag,std::complex<double> > LLSQLatticeMultiData; 


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


    std::vector<ADAT::Handle<LLSQLatticeMultiData> > build_multi_correlators(const ThreePointCorrIni_t &ini)
    {
      load(ini); 
      return build_multi_correlators(); 
    }

    std::vector<ADAT::Handle<LLSQLatticeMultiData> > build_multi_correlators(void);



    bool have_ini;
    ThreePointCorrIni_t m_ini;
  };

}








#endif /* BUILD_CORRELATORS_H */
