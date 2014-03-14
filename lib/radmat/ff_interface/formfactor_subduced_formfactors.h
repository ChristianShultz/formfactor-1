#ifndef FORMFACTOR_SUBDUCED_FORMFACTORS_H
#define FORMFACTOR_SUBDUCED_FORMFACTORS_H 

#include "formfactor.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_cubic_invariants.h"


namespace radmat
{
  struct SubducedFormFactorRecipe_t; 
  REGISTER_STRINGIFY_TYPE( SubducedFormFactorRecipe_t ); 

  struct SubducedFormFactorRecipe_t
    : public FormFactorRecipe_t
  {
    typedef HelicityFormFactorRecipe_t h_rep;   

    SubducedFormFactorRecipe_t(const h_rep &m, 
        const rHandle<CubicRep_p> &l, 
        const rHandle<CubicRep_p> &r,
        const std::string &left_table_id,
        const std::string &right_table_id)
      : hel(m) , lefty(l) , righty(r) , lt(left_table_id) , rt(right_table_id)
    { }

    SubducedFormFactorRecipe_t(const SubducedFormFactorRecipe_t &o)
      : hel(o.hel) , lefty(o.lefty) , righty(o.righty) , lt(o.lt) , rt(o.rt)
    { }

    SubducedFormFactorRecipe_t& operator=(const SubducedFormFactorRecipe_t &o)
    {
      if(this != &o)
      {
        hel = o.hel; 
        lefty = o.lefty; 
        righty = o.righty; 
        lt = o.lt;
        rt = o.rt; 
      }
      return *this;
    }

    virtual ~SubducedFormFactorRecipe_t() {}

    virtual int nFacs() const { return hel.nFacs(); }
    virtual std::string ff() const { return hel.id(); }
    std::map<int,std::string> ff_ids() const {return hel.ff_ids();}
    virtual rHandle<FFRep_p> left_rep() const { return FormFactorRecipe_t::call(lefty->rep_id());}
    virtual rHandle<FFRep_p> right_rep() const { return FormFactorRecipe_t::call(righty->rep_id());}
    virtual std::string id() const { return Stringify<SubducedFormFactorRecipe_t>(); }

    h_rep hel; 
    rHandle<CubicRep_p> lefty,righty ; 
    std::string lt,rt; 
  }; 




  struct SubducedFormFactor;
  REGISTER_STRINGIFY_TYPE(SubducedFormFactor); 

  struct SubducedFormFactor
    : public FormFactorBase_t
  {
    typedef FormFactorBase_t::recipe_h recipe_h; 

    SubducedFormFactor(const recipe_h &r)
      : FormFactorBase_t(r)
    { }

    SubducedFormFactor(const SubducedFormFactor &o)
      : FormFactorBase_t(o)
    { }

    SubducedFormFactor& operator=(const SubducedFormFactor &o)
    {
      FormFactorBase_t::operator=(o);
    }

    virtual ~SubducedFormFactor() {}

    virtual std::string id() const { return Stringify<SubducedFormFactor>(); }

    virtual itpp::Mat<std::complex<double> > 
      generate_ffs(const MomRowPair_t &lefty, 
          const MomRowPair_t &righty, 
          const double mom_fac) const; 

  }; 


  namespace SubducedFormFactorDecompositionFactoryEnv
  {
    bool registerAll(); 
    std::string build_id(const rHandle<CubicRep_p> &lefty,
        const rHandle<CubicRep_p> &righty); 
  }


} // radmat


#endif /* FORMFACTOR_SUBDUCED_FORMFACTORS_H */
