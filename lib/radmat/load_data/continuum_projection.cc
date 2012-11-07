/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : continuum_projection.cc

 * Purpose :

 * Creation Date : 24-10-2012

 * Last Modified : Wed Nov  7 10:39:41 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "continuum_projection.h"
#include "simple_world.h"
#include "invert_subduction.h"
#include "hadron/irrep_util.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "formfac/formfac_qsq.h" 
#include "radmat/utils/pow2assert.h"
#include "semble/semble_meta.h"
#include <string>

namespace radmat
{


  namespace 
  {
    Hadron::KeyHadronNPartNPtCorr_t::NPoint_t genNPt(const simpleWorld::ContinuumMatElem::State & cstate,
        const LatticeExprPrimitive & lstate)
    {
      Hadron::KeyHadronNPartNPtCorr_t::NPoint_t npt;
      npt.t_slice = cstate.t_slice;
      npt.irrep.row = lstate.row;
      npt.irrep.twoI_z = cstate.state.twoI_z;
      npt.irrep.mom = cstate.state.mom;
      npt.irrep.creation_op = cstate.state.creation_op;
      npt.irrep.smearedP = cstate.state.smearedP;
      npt.irrep.ops.resize(1);
      npt.irrep.ops[1].mom_type = FF::canonicalOrder(cstate.state.mom);
      npt.irrep.ops[1].name = cstate.state.name + std::string("_") + lstate.irrep;
      return npt; 
    }

  } // namespace anonomyous 



  // ---------------------------------------------------------------------------------



  ContinuumProjector::ContinuumProjector(const simpleWorld::ContinuumMatElem::State &obj)
    : m_init(false) 
  {
    load(obj);
  }

  void ContinuumProjector::load(const simpleWorld::ContinuumMatElem::State &state)
  {
    state_info = state;
    std::string LG = Hadron::generateLittleGroup(state_info.state.mom);
    ContinuumBosonExprPrimitive meson(state_info.state.J,
        state_info.state.parity,
        state_info.state.H,
        LG);
    ListLatticeIrrepExpr_t expr = invertSubduction(meson); 

    m_expr.m_expr.clear();

    ListLatticeIrrepExpr_t::const_iterator it;
    for(it = expr.begin(); it != expr.end(); ++it)
      m_expr.push_back(Expr_t(it->m_coeff, genNPt(state_info,it->m_obj)));

    m_init = true;
  } 






  // ---------------------------------------------------------------------------------




  namespace
  {
    Hadron::KeyHadronNPartNPtCorr_t genKey(const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t & sink, 
        const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t & ins,  
        const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t & source,
        const std::string & ensemble)
    {
      Hadron::KeyHadronNPartNPtCorr_t ret;
      ret.ensemble = ensemble; 
      ret.npoint.resize(3);
      ret.npoint[1] = sink;  // 1-based arrays are used here
      ret.npoint[2] = ins;
      ret.npoint[3] = source;

      return ret;
    }




  } // namespace anonomyous 



  ContinuumThreePointProjector::ContinuumThreePointProjector(const simpleWorld::ContinuumMatElem &elem)
  {
    load(elem);
  }


  void ContinuumThreePointProjector::load(const simpleWorld::ContinuumMatElem &elem)
  {
    POW2_ASSERT(elem.insertion.size() == 1); 
    ContinuumProjector source(elem.source);
    ContinuumProjector sink(elem.sink);
    ContinuumProjector ins(elem.insertion[0]);

    three_point_info = elem;

    m_expr.m_expr.clear();

    ContinuumProjector::const_iterator it_source, it_sink, it_ins;


    // this is basically summing over all cubic irreps and rows and weighting by non-zero subduction 
    // coefficients in order to "invert" the subduction procedure 
    for(it_source = source.begin(); it_source != source.end(); ++it_source)
      for(it_ins = ins.begin(); it_ins != ins.end(); ++it_ins)
        for(it_sink = sink.begin(); it_sink != sink.end(); ++it_sink)
          m_expr.push_back( Expr_t((it_source->m_coeff) * (it_ins->m_coeff) * (it_sink->m_coeff),
                genKey(it_source->m_obj, it_ins->m_obj, it_sink->m_obj, elem.ensemble)));

    m_init = true;
  }

} // namespace radmat

