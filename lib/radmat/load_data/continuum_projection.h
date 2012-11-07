#ifndef CONTINUUM_PROJECTION_H
#define CONTINUUM_PROJECTION_H


#include "simple_world.h"
#include "invert_subduction.h"
#include "hadron/hadron_npart_npt_corr.h"


namespace radmat
{

  struct ContinuumProjector
  {
    // usefull typedefs
    typedef ListLatticeIrrepExpr_t::Coeff_t Coeff_t;
    typedef Hadron::KeyHadronNPartNPtCorr_t::NPoint_t Obj_t;
    typedef ObjExpr_t<Coeff_t, Obj_t> Expr_t;
    typedef ListObjExpr_t<Coeff_t, Obj_t> List_t;
    typedef List_t::const_iterator const_iterator; 

    // constructors -- only simple world stuff for now but framework should be extensible
    ContinuumProjector(void) : m_init(false) {}
    ContinuumProjector(const simpleWorld::ContinuumMatElem::State &);

    // load from a simple world 
    void load(const simpleWorld::ContinuumMatElem::State &);

    // iteration for looping
    const_iterator begin(void) const {return m_expr.begin();}
    const_iterator end(void) const {return m_expr.end();}

    // peek
    const List_t& expr(void) const {return m_expr;}


    private:
    List_t m_expr;
    simpleWorld::ContinuumMatElem::State state_info; 
    bool m_init;
  };



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


  struct ContinuumThreePointProjector
  {
    // usefull typedefs
    typedef ContinuumProjector::Coeff_t Coeff_t;
    typedef Hadron::KeyHadronNPartNPtCorr_t Obj_t;
    typedef ObjExpr_t<Coeff_t, Obj_t> Expr_t;
    typedef ListObjExpr_t<Coeff_t,Obj_t> List_t;
    typedef List_t::const_iterator const_iterator;


    // constructors -- only simple world stuff for now but framework should be extensible
    ContinuumThreePointProjector(void) : m_init(false) {}
    ContinuumThreePointProjector(const simpleWorld::ContinuumMatElem &);

    // load from a simple world
    void load(const simpleWorld::ContinuumMatElem &);

    // iteration for looping
    const_iterator begin(void) const {return m_expr.begin();}
    const_iterator end(void) const {return m_expr.end();}

    // peek
    const List_t & expr(void) const {return m_expr;}


    private:

    List_t m_expr;
    simpleWorld::ContinuumMatElem three_point_info;
    bool m_init;  
  };


} // namespace radmat

















#endif /* CONTINUUM_PROJECTION_H */
