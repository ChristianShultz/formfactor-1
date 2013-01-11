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
    redstarCircularMatElem_t(const simpleWorld::ContinuumMatElem &, const std::string &source_id, const std::string &sink_id); 

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


  // the lattice matrix elements are in a helicity basis but the linear system
  // solvers are set up to use a cartesian basis.. use this to do a trivial change of basis
  struct redstarCartesianMatElem_p
  {
    // save some typing and in the process make this virtually unreadable to anyone other than myself..
    typedef redstarCircularMatElem_t::redKey redKey;
    typedef redstarCircularMatElem_t::normKey normKey;
    typedef redstarCircularMatElem_t::subducedOp subducedOp;
    typedef redstarCircularMatElem_t::listSubducedOp listSubducedOp;
    typedef redstarCircularMatElem_t::subducedInsertion subducedInsertion;
    typedef redstarCircularMatElem_t::listSubducedInsertion listSubducedInsertion;


    redstarCartesianMatElem_p(const redstarCircularMatElem_t &elem) 
    {
      makeCartesian(elem);
    }

    redstarCartesianMatElem_p(const simpleWorld::ContinuumMatElem &a, const std::string &source_id, const std::string &sink_id)
    {
      redstarCircularMatElem_t dum(a,source_id,sink_id);
      ensemble = a.ensemble; 
      makeCartesian(dum);
    }

    // !! NB:
    // project out the x,y components from the +,- helicity components .. assumed the usual normalization but 
    // this should be checked for consistency in case theres some strange lattice thingy I don't know about 
    void makeCartesian(const redstarCircularMatElem_t &tmp)
    {
      m_source = tmp.m_source;
      m_sink = tmp.m_sink;

      m_t = tmp.m_time;
      m_z = tmp.m_zero;
      m_x = std::pair<bool,listSubducedInsertion>(tmp.m_plus.first && tmp.m_minus.first,
          SEMBLE::toScalar(std::complex<double>(1./sqrt(2.),0.))*(tmp.m_plus.second + tmp.m_minus.second));
      m_y = std::pair<bool,listSubducedInsertion>(tmp.m_plus.first && tmp.m_minus.first,
          SEMBLE::toScalar(std::complex<double>(0.,-1./sqrt(2.)))*(tmp.m_plus.second - tmp.m_minus.second));

    }

    listSubducedOp m_source;
    listSubducedOp m_sink;
    std::pair<bool,listSubducedInsertion> m_t; 
    std::pair<bool,listSubducedInsertion> m_x;
    std::pair<bool,listSubducedInsertion> m_y;
    std::pair<bool,listSubducedInsertion> m_z;
    std::string ensemble;  
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

    redstarCartMatElemLorentzComponent get_component(const int lorentz_index)
    {
      POW2_ASSERT(lorentz_components.find(lorentz_index) != lorentz_components.end());
      return lorentz_components.find(lorentz_index)->second; 
    }

    std::map<int,redstarCartMatElemLorentzComponent> lorentz_components; 

  };



} // namespace radmat 



#endif /* GENERATE_REDSTAR_XML_H */
