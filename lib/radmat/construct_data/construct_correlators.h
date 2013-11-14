#ifndef CONSTRUCT_CORRELATORS_H
#define CONSTRUCT_CORRELATORS_H

/*
 * This class is responsible for taking the ini and generating either
 * the redstar xml we need to generate or for combining lattice data 
 * in order to produce correlation functions corresponding to matrix
 * elements in a cartesian basis.
 */

#include "construct_correlators_xml.h" 
#include "invert_subduction.h"
#include "lattice_multi_data_object.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "adat/handle.h"
#include <vector>

namespace radmat
{
  struct ConstructCorrelators
  {
    // some boiler plate stuff
    ConstructCorrelators(void) : have_ini(false) {}
    ConstructCorrelators(const ThreePointCorrIni_t &ini) : have_ini(true) , m_ini(ini) {}

    // load an ini
    void load(const ThreePointCorrIni_t &ini)
    {
      registerSubductionTables(); 
      have_ini = true;
      m_ini = ini; 

      print_redstar(m_ini.threePointCorrXMLIni.redstar); 
    }

    void print_redstar(const AbstractMergeNamedObject &); 

    // sum over irreps to produce cont correlators
    std::vector<ADAT::Handle<LLSQLatticeMultiData> >
     construct_multi_correlators(const ThreePointCorrIni_t &ini)
    {
      load(ini); 
      return construct_multi_correlators(); 
    }

    std::vector<ADAT::Handle<LLSQLatticeMultiData> > 
      construct_multi_correlators(void) const;

    // figure out what xml we would need
    std::vector<Hadron::KeyHadronNPartNPtCorr_t> 
      construct_correlator_xml(const ThreePointCorrIni_t &ini)
    {
      load(ini);
      return construct_correlator_xml(); 
    }

    std::vector<Hadron::KeyHadronNPartNPtCorr_t>  
      construct_correlator_xml(void) const;
  
    // public data store
    bool have_ini;
    ThreePointCorrIni_t m_ini;
  };

} // radmat



#endif /* CONSTRUCT_CORRELATORS_H */
