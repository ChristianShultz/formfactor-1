#ifndef GENERATE_REDSTAR_XML_H
#define GENERATE_REDSTAR_XML_H


#include "radmat/load_data/simple_world.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/pow2assert.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "radmat_overlap_key_val_db.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"

#include <vector>
#include <utility>

namespace radmat
{




  // struct to hold the xml and coefficients associated with a given lattice xml mat elem
  struct redstarCircularMatElem_t
  {
    typedef Hadron::KeyHadronNPartNPtCorr_t::NPoint_t redKey; 
    typedef RadmatExtendedKeyHadronNPartIrrep_t normKey;

    struct OperatorKeyData
    {
      OperatorKeyData(const redKey &red, const normKey &norm)
        : operatorKey(red) , normalizationKey(norm) {}

      redKey operatorKey;
      normKey normalizationKey; 
    };


    typedef ObjExpr_t<ENSEM::Complex,OperatorKeyData> subducedOp; 
    typedef ListObjExpr_t<ENSEM::Complex,OperatorKeyData> listSubducedOp; 

    typedef ObjExpr_t<ENSEM::Complex,redKey> subducedInsertion;
    typedef ListObjExpr_t<ENSEM::Complex,redKey> listSubducedInsertion;


    // overload this constructor to deal with optimized operators???
    redstarCircularMatElem_t(const simpleWorld::ContinuumMatElem &,
        const std::string &source_id,
        const std::string &sink_id); 

    // the extra name of the state for the normalization database
    const std::string m_source_id;
    const std::string m_sink_id; 

    // lists of subduced operators and associated coefficients 
    listSubducedOp m_source;
    listSubducedOp m_sink;
    std::pair<bool,listSubducedInsertion> m_time; 
    std::pair<bool,listSubducedInsertion> m_plus;
    std::pair<bool,listSubducedInsertion> m_zero;
    std::pair<bool,listSubducedInsertion> m_minus; 
  };


  // the actual real world interface 
  struct redstarCartMatElem
  {

    // save some typing and in the process make this virtually unreadable to anyone other than myself..
    typedef redstarCircularMatElem_t::redKey redKey;
    typedef redstarCircularMatElem_t::normKey normKey;
    typedef redstarCircularMatElem_t::subducedOp subducedOp;
    typedef redstarCircularMatElem_t::listSubducedOp listSubducedOp;
    typedef redstarCircularMatElem_t::subducedInsertion subducedInsertion;
    typedef redstarCircularMatElem_t::listSubducedInsertion listSubducedInsertion;

    // this is the real world guy that we will use when querrying the database
    struct redstarCartMatElemLorentzComponent
    {

      // hold the key for the three point corr and the two keys associated with energy/overlap on the source and sink
      struct threePointKey
      {

        threePointKey(const Hadron::KeyHadronNPartNPtCorr_t &red_xml, 
            const RadmatExtendedKeyHadronNPartIrrep_t  &source_xml,
            const RadmatExtendedKeyHadronNPartIrrep_t &sink_xml)
          : redstar_xml(red_xml), source_normalization(source_xml) , sink_normalization(sink_xml)  { }

        Hadron::KeyHadronNPartNPtCorr_t redstar_xml;
        RadmatExtendedKeyHadronNPartIrrep_t source_normalization;
        RadmatExtendedKeyHadronNPartIrrep_t sink_normalization;
      };

      // make a list of these keys with a given weight determined by the subduction procedure 
      typedef ObjExpr_t<ENSEM::Complex,threePointKey> threePointXMLKey;
      typedef ListObjExpr_t<ENSEM::Complex,threePointKey> listThreePointXMLKey; 
      typedef listThreePointXMLKey::const_iterator const_iterator;  

      redstarCartMatElemLorentzComponent(void) : active(false) {  }

      // even though this is quantum mechanics we are going to read from left to right..
      redstarCartMatElemLorentzComponent(const listSubducedOp &source,
          const std::pair<bool,listSubducedInsertion> &ins,
          const listSubducedOp &sink,
          const std::string &ensemble);

      // stl like overloads for looping once we start spaming database fetches 
      const_iterator begin(void) const {return m_list.begin();}
      const_iterator end(void) const {return m_list.end();}
      bool is_active(void) const {return active;}

      bool active; 
      listThreePointXMLKey m_list; 
    };  // redstarCartMatElemLorentzComponent


    // make a list of these keys with a given weight determined by the subduction procedure 
    typedef redstarCartMatElemLorentzComponent::threePointXMLKey threePointXMLKey;
    typedef redstarCartMatElemLorentzComponent::listThreePointXMLKey listThreePointXMLKey; 
    typedef listThreePointXMLKey::const_iterator const_iterator;


    redstarCartMatElem(const simpleWorld::ContinuumMatElem & , const std::string &source_id, const std::string &sink_id);

    ~redstarCartMatElem(void); 

    redstarCartMatElemLorentzComponent get_component(const int lorentz_index)
    {
      POW2_ASSERT(lorentz_components.find(lorentz_index) != lorentz_components.end());
      return lorentz_components.find(lorentz_index)->second; 
    }


    simpleWorld::ContinuumMatElem m_elem; 
    std::map<int,redstarCartMatElemLorentzComponent> lorentz_components; 

  };



} // namespace radmat 



#endif /* GENERATE_REDSTAR_XML_H */
