#ifndef BUILD_CORRELATORS_H
#define BUILD_CORRELATORS_H

/*
 * This class is responsible for taking the ini and generating either
 * the redstar xml we need to generate or for combining lattice data 
 * in order to produce correlation functions corresponding to matrix
 * elements in a cartesian basis.
 */

#include "g_parity_world_xml.h"
#include "g_parity_world.h"
#include "g_parity_world_generate_redstar_xml.h"
#include "radmat_database_interface.h"
#include "build_correlators_xml.h" 
#include "lattice_multi_data_tag.h"
#include "lattice_multi_data_object.h"
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
