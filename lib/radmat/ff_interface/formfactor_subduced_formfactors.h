#ifndef FORMFACTOR_SUBDUCED_FORMFACTORS_H
#define FORMFACTOR_SUBDUCED_FORMFACTORS_H 

#include "formfactor.h"
#include "radmat/utils/obj_expr_t.h"
#include "ensem/ensem.h"


namespace radmat
{
  struct SubducedFormFactorRecipe_t; 
  REGISTER_STRINGIFY_TYPE( SubducedFormFactorRecipe_t ); 

  struct SubducedFormFactorRecipe_t
    : public FormFactorRecipe_t
  {
    typedef rHandle<LorentzFFAbsBase_t> mat_h;

    SubducedFormFactorRecipe_t(const mat_h &m, const rHandle<CubicRep_p> &l, const rHandle<CubicRep_p> &r)
      : mat(m) , lefty(l) , righty(r)
    { }

    SubducedFormFactorRecipe_t(const SubducedFormFactorRecipe_t &o)
      : mat(o.mat) , lefty(o.lefty) , righty(o.righty)
    { }

    SubducedFormFactorRecipe_t& operator=(const SubducedFormFactorRecipe_t &o)
    {
      if(this != &o)
      {
        mat = o.mat; 
        lefty = o.lefty; 
        righty = o.righty; 
      }
      return *this;
    }

    virtual ~SubducedFormFactorRecipe_t() {}

    virtual int nFacs() const { return mat->nFacs(); }
    virtual std::string ff() const { return mat->id(); }
    std::map<int,std::string> ff_ids() const {return mat->ff_ids();}
    virtual int right_spin() const { return mat->right_spin(); }
    virtual int left_spin() const { return mat->left_spin(); }
    virtual int right_row() const { return righty->rep_row(); }
    virtual int left_row() const { return lefty->rep_row(); }

    mat_h mat; 
    rHandle<CubicRep_p> lefty,righty ; 
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
