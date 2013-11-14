#ifndef CONSTRUCT_CORRELATORS_UTILS_H
#define CONSTRUCT_CORRELATORS_UTILS_H 


#include "construct_correlators_xml.h"
#include "lattice_multi_data_tag_redstar_interface.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "radmat_database_interface.h"
#include "radmat_overlap_key_val_db.h"
#include "ensem/ensem.h"
#include "adat/map_obj.h"
#include <vector>
#include <string>
#include <map> 
#include <utility>

namespace radmat
{

  std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > 
    sort_tagged_corrs_by_Q2(const std::vector<TaggedEnsemRedstarNPtBlock> &); 


  // the database type we will be using 
  typedef radmatAllConfDatabaseInterface< Hadron::KeyHadronNPartNPtCorr_t,
          ENSEM::EnsemVectorComplex,
          RadmatExtendedKeyHadronNPartIrrep_t,
          RadmatMassOverlapData_t> DatabaseInterface_t;

  // NB hardwire of types from above in here
  struct DatabaseInterface_k
  {    
    DatabaseInterface_k(const Hadron::KeyHadronNPartNPtCorr_t &k ,
        const std::string &id_sink, const std::string &id_source)
      : npt(k)
    {
      sink = RadmatExtendedKeyHadronNPartIrrep_t(id_sink,k.npoint[1].irrep);
      source = RadmatExtendedKeyHadronNPartIrrep_t(id_source,k.npoint[3].irrep); 
    }

    Hadron::KeyHadronNPartNPtCorr_t npt;
    RadmatExtendedKeyHadronNPartIrrep_t sink,source;  
  };  


  struct ConstructCorrsMatrixElement
  {
    ConstructCorrsMatrixElement(const ENSEM::EnsemVectorComplex &d, 
        const LatticeMultiDataTag &t,
        const bool s
        )
      : data(d) , tag(t) , success(s)
    {  }

    ENSEM::EnsemVectorComplex data; 
    LatticeMultiDataTag tag; 
    bool success; 
  };

  std::pair<bool,std::vector<ConstructCorrsMatrixElement> >
    build_correlators(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        std::vector<Hadron::KeyHadronNPartNPtCorr_t> &missed_xml, 
        std::vector<RadmatExtendedKeyHadronNPartIrrep_t> &missed_norm, 
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &,
        const DatabaseInterface_t & );


} // radmat


#endif /* CONSTRUCT_CORRELATORS_UTILS_H */