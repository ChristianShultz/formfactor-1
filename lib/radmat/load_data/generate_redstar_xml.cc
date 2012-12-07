/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : generate_redstar_xml.cc

 * Purpose :

 * Creation Date : 03-12-2012

 * Last Modified : Tue Dec  4 13:45:15 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "generate_redstar_xml.h"
#include "radmat/load_data/invert_subduction.h"
#include "radmat/load_data/simple_world.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/pow2assert.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "hadron/hadron_npart_irrep.h"
#include "hadron/irrep_util.h"
#include "semble/semble_key_val_db.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"; 

#include <sstream>

#if 0



// struct to hold the xml and coefficients associated with a given lattice xml mat elem
struct redstarCircularMatElem_t
{
  typedef Hadron::KeyHadronNPartNPtCorr_t::NPoint_t redKey; 
  typedef SEMBLE::SembleExtendedKeyHadronNPartIrrep_t normKey;

  struct OperatorKeyData
  {
    OperatorKeyData(const redKey &red, const normKey &norm)
      : operatorKey(red) , normalizationKey(norm) {}

    redKey operatorKey;
    normKey normalizationKey; 
  };


  typedef ObjExpr_t<ENSEM::Complex,OperatorKeyData> subducedOp; 
  typedef ListObjExpr_t<ENSEM::Complex,subducedOp> listSubducedOp; 

  typedef ObjExpr_t<ENSEM::Complex,redKey> subducedInsertion;
  typedef ListObjExpr_t<ENSEM::Complex,subducedInsertion> listSubducedInsertion;


  // overload this constructor to deal with optimized operators???
  redstarCircularMatElem_t(const simpleWorld::ContinuumMatElem &, const std::string &source_id, const std::string &sink_id); 

  // the extra name of the state for the normalization database
  const std::string m_source_id;
  const std::string m_sink_id; 

  // lists of subduced operators and associated coefficients 
  listSubducedOp m_source;
  listSubducedOp m_sink
    std::pair<bool,listSubducedInsertion> m_time; 
  std::pair<bool,listSubducedInsertion> m_plus;
  std::pair<bool,listSubducedInsertion> m_zero;
  std::pair<bool,listSubducedInsertion> m_minus; 
};
#endif


// CIRCULAR BASIS STUFF
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

namespace radmat
{

  namespace 
  {
    // save some typing and in the process make this virtually unreadable to anyone other than myself..
    typedef redstarCircularMatElem_t::redKey redKey;
    typedef redstarCircularMatElem_t::normKey normKey;
    typedef redstarCircularMatElem_t::subducedOp subducedOp;
    typedef redstarCircularMatElem_t::listSubducedOp listSubducedOp;
    typedef redstarCircularMatElem_t::subducedInsertion subducedInsertion;
    typedef redstarCircularMatElem_t::listSubducedInsertion listSubducedInsertion;


    listSubducedOp getLatticeSubducedOp(const simpleWorld::ContinuumStatePrimitive &cont,
        const std::string &pid,
        const int t_slice)
    {
      listSubducedOp ret; 
      ContinuumBosonExprPrimitive boson(cont.J,cont.parity,cont.H,Hadron::generateLittleGroup(cont.mom)); 
      ListLatticeIrrepExpr_t lattice_boson = invertSubduction(boson); 
      ListLatticeIrrepExpr_t::const_iterator it;

      for(it = lattice_boson.begin(); it != lattice_boson.end(); ++it)
      {
        std::stringstream concat_name; 
        concat_name << cont.name;

        if(boson.group == "Oh")
        {
          concat_name << "_" << it->m_obj.irrep;
        }
        else
        {
          concat_name << "_H" << boson.H << it->m_obj.group << it->m_obj.irrep;
        }

        Hadron::KeyHadronNPartIrrep_t base; 
        base.ops.resize(1); 
        base.ops[1].name = concat_name.str(); 
        base.ops[1].mom_type = FF::canonicalOrder(cont.mom); 
        base.row = it->m_obj.row; 
        base.twoI_z = cont.twoI_z; 
        base.mom = cont.mom;
        base.creation_op = cont.creation_op;
        base.smearedP = cont.smearedP;

        Hadron::KeyHadronNPartNPtCorr_t::NPoint_t npt;
        npt.t_slice = t_slice;
        npt.irrep = base;   

        SEMBLE::SembleExtendedKeyHadronNPartIrrep_t norm(pid,base);     
        redstarCircularMatElem_t::OperatorKeyData key(npt,norm);

        ret.push_back(subducedOp(it->m_coeff,key));
      }   
      return ret; 
    }


    listSubducedOp getLatticeSubducedOp(const simpleWorld::ContinuumMatElem::State &state, const std::string &pid)
    {
      return getLatticeSubducedOp(state.state,pid,state.t_slice);
    }


    struct circLorentzInsertion
    {
      circLorentzInsertion(void)
        : time(std::pair<bool,listSubducedInsertion>(false,listSubducedInsertion())) , plus(time), zero(time) , minus(time) {}

      std::pair<bool,listSubducedInsertion> time;
      std::pair<bool,listSubducedInsertion> plus;
      std::pair<bool,listSubducedInsertion> zero;
      std::pair<bool,listSubducedInsertion> minus; 
    };

    circLorentzInsertion getLatticeSubducedInsertion(const simpleWorld::ContinuumInsertion &ins, const std::string &pid)
    {
      circLorentzInsertion ret;
      listSubducedOp work;
      listSubducedOp::const_iterator work_it;
      std::map<std::string,simpleWorld::ContinuumInsertion::op_insertion>::const_iterator it;

      for(it = ins.insertion_map.begin(); it != ins.insertion_map.end(); ++it)
      {
        work = getLatticeSubducedOp(it->second,pid,ins.t_slice);

        listSubducedInsertion ins_list;

        for(work_it = work.begin(); work_it != work.end(); ++work_it)
          ins_list.push_back(subducedInsertion(work_it->m_coeff,work_it->m_obj.operatorKey));

        if(it->first == "t")
          ret.time = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else if(it->first == "p")
          ret.plus = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else if(it->first == "0")
          ret.zero = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else if(it->first == "m")
          ret.minus = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else
        {
          std::cout << __func__ << ": error: can't match to a circular lorentz index" << std::endl;
          exit(1);
        }

      }

      return ret;
    }


  } // namespace anonomyous 



  redstarCircularMatElem_t::redstarCircularMatElem_t(const simpleWorld::ContinuumMatElem & cont,
      const std::string &source_id, 
      const std::string &sink_id)
    : m_source_id(source_id), m_sink_id(sink_id)
  {
    m_source = getLatticeSubducedOp(cont.source,m_source_id);
    m_sink = getLatticeSubducedOp(cont.sink,m_sink_id);
    circLorentzInsertion foo = getLatticeSubducedInsertion(cont.insertion,std::string("foobar"));
    m_time = foo.time;
    m_plus = foo.plus;
    m_zero = foo.zero;
    m_minus = foo.minus;
  }

} // namespace radmat




// CARTESIAN BASIS STUFF
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

namespace radmat
{

  namespace
  {
    redstarCartMatElem::redstarCartMatElemLorentzComponent::threePointKey genThreePointKey(const listSubducedOp::ListObj_t &source, 
        const listSubducedInsertion::ListObj_t &ins,
        const listSubducedOp::ListObj_t &sink)
    {
      Hadron::KeyHadronNPartNPtCorr_t redstar_xml;
      redstar_xml.npoint.resize(3); // 1 based array
      redstar_xml.npoint[1] = source.m_obj.operatorKey;
      redstar_xml.npoint[2] = ins.m_obj;
      redstar_xml.npoint[3] = sink.m_obj.operatorKey;

      return redstarCartMatElem::redstarCartMatElemLorentzComponent::threePointKey(redstar_xml,source.m_obj.normalizationKey,sink.m_obj.normalizationKey); 
    }

  } // namespace anonomyous


  redstarCartMatElem::redstarCartMatElemLorentzComponent::redstarCartMatElemLorentzComponent(const listSubducedOp &source,
        const std::pair<bool,listSubducedInsertion> &ins,
        const listSubducedOp &sink) : active(false)
  {
    // leave the list empty if the init variable from the xml was false for whatever reason
    if(ins.first)
    {
      listSubducedOp::const_iterator it_source,it_sink;
      listSubducedInsertion::const_iterator it_ins;

      // run the triple for loop and keep track of the subduction coefficients 
      for(it_source = source.begin(); it_source != source.end(); ++it_source)
        for(it_ins = ins.second.begin(); it_ins != ins.second.end(); ++it_ins)
          for(it_sink = sink.begin(); it_sink != sink.end(); ++it_sink)
            m_list.push_back(threePointXMLKey(it_source->m_coeff * it_ins->m_coeff * it_sink->m_coeff,
                  genThreePointKey(*it_source,*it_ins,*it_sink)
                  )
                );
      active = true; 
    }

  }


  redstarCartMatElem::redstarCartMatElem(const simpleWorld::ContinuumMatElem &cont,
       const std::string &source_id,
       const std::string &sink_id)
  {
    redstarCartesianMatElem_p foobar(cont,source_id,sink_id);

    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(0,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_t,foobar.m_sink))
        );
    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(1,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_x,foobar.m_sink))
        );
    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(2,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_y,foobar.m_sink))
        );
    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(3,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_z,foobar.m_sink))
        );

    if(foobar.m_t.first)
      lorentz_components[0].active = true;
    if(foobar.m_x.first)
      lorentz_components[1].active = true;
    if(foobar.m_y.first)
      lorentz_components[2].active = true;
    if(foobar.m_z.first)
      lorentz_components[3].active = true;
    
  }


} // namespace radmat

