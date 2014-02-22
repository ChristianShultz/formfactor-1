#ifndef FORMFACTOR_HELICITY_FORMFACTORS_H
#define FORMFACTOR_HELICITY_FORMFACTORS_H 




#include "formfactor.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_factory.h"
#include "formfactor_spherical_invariants.h"
#include "radmat/ff/lorentzff_formfactor_abs_base_cfg.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"



namespace radmat
{

  struct HelicityFormFactorRecipe_t;
  REGISTER_STRINGIFY_TYPE( HelicityFormFactorRecipe_t ); 

  struct HelicityFormFactorRecipe_t
  {
    typedef rHandle<LorentzFFAbsBase_t> mat_h; 
    typedef rHandle<SpherRep_p> rep_h; 

    HelicityFormFactorRecipe_t(mat_h &f, rep_h &l, rep_h &r)
      : mat(f) , lefty(l) , righty(r)
    {}

    HelicityFormFactorRecipe_t( const HelicityFormFactorRecipe_t &o)
      : mat(o.mat) , lefty(o.lefty) , righty(o.righty)
    {}

    HelicityFormFactorRecipe_t& operator=( const HelicityFormFactorRecipe_t &o)
    {
      if (this != &o) 
      {
        mat = o.mat; 
        lefty = o.lefty;
        righty = o.righty; 
      }
      return *this;
    }

    virtual ~HelicityFormFactorRecipe_t() {}

    virtual int nFacs() const { return mat->nFacs(); }
    virtual std::string ff() const { return mat->id(); }
    virtual std::map<int,std::string> ff_ids() const { return mat->ff_ids(); }
    virtual int right_spin() const { return mat->right_spin(); }
    virtual int left_spin() const { return mat->left_spin(); }
    virtual int right_row() const { return righty->rep_row(); }
    virtual int left_row() const { return lefty->rep_row(); }
    virtual std::string id() const { return Stringify<HelicityFormFactorRecipe_t>(); }

    rHandle<LorentzFFAbsBase_t> mat; 
    rHandle<SpherRep_p> lefty, righty; 
  };


  struct HelicityFormFactor; 
  REGISTER_STRINGIFY_TYPE( HelicityFormFactor ); 

  struct HelicityFormFactor
    : public FormFactorBase_t
  {
    typedef FormFactorBase_t::recipe_h recipe_h;  

    HelicityFormFactor( const recipe_h & recipe_)
      : FormFactorBase_t(recipe)
    { }

    HelicityFormFactor& operator=(const HelicityFormFactor &o)
    {
      FormFactorBase_t::operator=(o); 
    }

    HelicityFormFactor( const HelicityFormFactor &o)
      : FormFactorBase_t(o)
    { }

    virtual ~HelicityFormFactor() {}

    virtual itpp::Mat<std::complex<double> >
      generate_ffs( const MomRowPair_t &lefty, 
         const MomRowPair_t &righty, 
         const double mom_fac) const; 
  };


  namespace HelicityFormFactorDecompositionFactoryEnv
  {
    bool registerAll(void) ;
    std::string build_id( const rHandle<SpherRep_p> &lefty,
        const rHandle<SpherRep_p> &righty ); 
  }


} // radmat

#endif /* FORMFACTOR_HELICITY_FORMFACTORS_H */
