#ifndef BUILD_CORRELATORS_H
#define BUILD_CORRELATORS_H


#include "g_parity_world_xml.h"
#include "g_parity_world.h"
#include "g_parity_world_generate_redstar_xml.h"
#include "radmat_database_interface.h"
#include "build_correlators_xml.h" 
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


  // this is a tag for a row in the linear system 
  struct LatticeMultiDataTag
  {
    //! Constructor
    LatticeMultiDataTag(void);

    LatticeMultiDataTag& operator=(const LatticeMultiDataTag &o); 

    //! ensemble qsquared
    ENSEM::EnsemReal Q2(void) const;

    //! splash
    void print_me(void) const;

    //! got sick of typing this
    std::string mom_string(void) const;

    //! energies
    std::string E_string(void) const;

    //! the value of Q2 we use when sorting and labeling
    void set_qsq_label(const double &q2) {qsq_label = q2;}

    //! the value of Q2 we use when sorting and labeling 
    double get_qsq_label(void) const {return qsq_label;}

    std::string splash_tag(void) const;

    // tags

    std::string file_id; // some unique string telling us what this is

    // for llsq system
    double qsq_label;        
    int jmu;                 
    std::string mat_elem_id; 
    Array<int> p_f;
    Array<int> p_i;
    ENSEM::EnsemReal E_f;
    ENSEM::EnsemReal E_i;
    double mom_fac; 
  };


  // the definition of a linear system that we are creating 
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


    std::vector<ADAT::Handle<LLSQLatticeMultiData> >
     build_multi_correlators(const ThreePointCorrIni_t &ini)
    {
      load(ini); 
      return build_multi_correlators(); 
    }

    std::vector<ADAT::Handle<LLSQLatticeMultiData> > build_multi_correlators(void);


    std::vector<Hadron::KeyHadronNPartNPtCorr_t> 
      build_correlator_xml(const ThreePointCorrIni_t &ini)
    {
      load(ini);
      return build_correlator_xml(); 
    }

    std::vector<Hadron::KeyHadronNPartNPtCorr_t>  build_correlator_xml(void);
  
    bool have_ini;
    ThreePointCorrIni_t m_ini;
  };

} // radmat



#endif /* BUILD_CORRELATORS_H */
