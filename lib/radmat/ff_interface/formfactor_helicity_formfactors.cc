/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_helicity_formfactors.cc

 * Purpose :

 * Creation Date : 22-02-2014

 * Last Modified : Wed 26 Feb 2014 01:39:28 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "formfactor.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_factory.h"
#include "formfactor_spherical_invariants.h"
#include "radmat/ff/lorentzff_canonical_cont_spin_formfactors.h"
#include "radmat/ff/lorentzff_formfactor_factory.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include <sstream>
#include <algorithm>

#define CHECK_HELICITY_EXPLICIT

typedef radmat::TheFormFactorInjectionRecipeCookbook Cookbook; 
typedef radmat::TheFormFactorInjectionRecipeFactory Factory; 

namespace radmat
{

  namespace
  {
    // sanity
#ifdef CHECK_HELICITY_EXPLICIT
    bool check_helicity(const int h, const int row, const int J)
    {
      return h == ( J - row + 1 ); 
    }
#endif

    struct print_elem_reg
    {
      static void print(const std::string &msg)
      {std::cout << msg << std::endl;}
    };

    struct dont_print_elem_reg
    {
      static void print(const std::string &msg) {}
    };
  }

  // use whatever is in the lorentzff_canonical_cont_spin decomp
  itpp::Mat<std::complex<double> >
    HelicityFormFactor::generate_ffs( const MomRowPair_t &lefty, 
        const MomRowPair_t &righty, 
        const double mom_fac) const
    {
      FormFactorBase_t::recipe_h recipe = FormFactorBase_t::get_recipe(); 
      POW2_ASSERT( recipe->id() == Stringify<HelicityFormFactorRecipe_t>() ); 
      const HelicityFormFactorRecipe_t * helicity_recipe;
      helicity_recipe = dynamic_cast< const HelicityFormFactorRecipe_t* >( recipe.get_ptr() );

#ifdef CHECK_HELICITY_EXPLICIT
      check_helicity( lefty.second, helicity_recipe->left_row() , helicity_recipe->left_spin() ); 
      check_helicity( righty.second, helicity_recipe->right_row(), helicity_recipe->right_spin() ); 
#endif

      return helicity_recipe->mat->operator()(lefty,righty,mom_fac); 
    }


  // now worry about dumping recipes into the factory
  namespace HelicityFormFactorDecompositionFactoryEnv
  {

    // local utility
    namespace
    { 
      // build the id from the lorentz spin factory, 
      // this is a class name like J0pJ3m
      std::string build_lorentz_spin_id(const rHandle<SpherRep_p> &lefty, 
          const rHandle<SpherRep_p> &righty)
      {
        std::stringstream ss; 
        ss << lefty->rep_id() << righty->rep_id();
        return ss.str(); 
      }

      // snarf (or slurp?) all continuum reps
      std::vector<rHandle<SpherRep_p> > 
        get_all_cont_reps()
        {
          std::vector<std::string> keys;
          std::vector<std::string>::const_iterator it; 
          std::vector<rHandle<SpherRep_p> > vals; 

          keys = ::radmat::SpherInvariantsFactoryEnv::all_keys(); 
          for(it = keys.begin(); it != keys.end(); ++it)
            vals.push_back( 
                ::radmat::SpherInvariantsFactoryEnv::callFactory(*it) ); 

          return vals;
        }

      // reg function
      typedef std::pair<rHandle<SpherRep_p> , rHandle<SpherRep_p> > rep_pair;

      FormFactorRecipe_t* gen_recipe( const std::string &s, const rep_pair &reps )
      {
        rHandle<LorentzFFAbsBase_t> mat; 
        mat = radmat::LorentzffFormFactorDecompositionFactoryEnv::callFactory(s); 
        return new HelicityFormFactorRecipe_t(mat,reps.first,reps.second); 
      }

      FFAbsBase_t* callback( const std::string &recipe_id )
      {
        PtrRecipeHolder::map_t::iterator r;  

        r = Cookbook::Instance().mappy.find(recipe_id); 
        bool success = true; 

        if( r == Cookbook::Instance().mappy.end() )
        {
          printer_function<console_print>( "missed " + recipe_id); 
          throw std::string("recipe_error"); 
        }

        if( r->second->id() != radmat::Stringify<HelicityFormFactorRecipe_t>())
        {
          printer_function<console_print>( "expected a " + Stringify<HelicityFormFactorRecipe_t>());
          printer_function<console_print>( "got a " + r->second->id());
          throw std::string("recipe_error"); 
        }

        HelicityFormFactorRecipe_t *ptr = dynamic_cast<HelicityFormFactorRecipe_t*>( r->second ); 
        rHandle<FormFactorRecipe_t> recipe( new HelicityFormFactorRecipe_t(*ptr) ); 

        return new HelicityFormFactor( recipe ); 
      }

      bool do_reg( const std::string &s, const rep_pair &r)
      {
        // the actual string is between the <> eg  <reg_id>
        printer_function<dont_print_elem_reg>(std::string("regged <") + s + std::string(">"));
        Cookbook::Instance().mappy.insert(std::make_pair(s, gen_recipe(s,r)) ); 
        bool b = true; 
        b = Factory::Instance().registerObject(s,callback);  
        if( !!! b ) 
          printer_function<console_print>( s + " failed to register" ); 
        return b; 
      }

    } // anonomyous



    // the factory reg id for helicity form factors
    std::string build_id( const rHandle<SpherRep_p> &lefty, 
        const rHandle<SpherRep_p> &righty)
    {
      std::stringstream ss; 
      ss << lefty->reg_id() << "__" << righty->reg_id();
      return ss.str(); 
    }

    // local registration flag
    namespace 
    {
      bool registered = false; 
    }

    // pump the factory 
    bool registerAll()
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if ( !!! registered )
      {
        typedef std::pair<rHandle<SpherRep_p> , rHandle<SpherRep_p> > rep_pair;
        typedef std::map<std::string, rep_pair > map_t; 
        typedef map_t::value_type value_type; 

        map_t ff_spin_keys;
        map_t::const_iterator it; 
        std::vector<std::string> ff_allowed_spin_keys; 
        std::vector<rHandle<SpherRep_p> > cont_rep_keys; 
        std::vector<rHandle<SpherRep_p> >::const_iterator i,j; 

        cont_rep_keys = get_all_cont_reps();

        for( i = cont_rep_keys.begin(); i != cont_rep_keys.end(); ++i)
          for( j = cont_rep_keys.begin(); j != cont_rep_keys.end(); ++j)
          {
            std::string s = build_lorentz_spin_id(*i,*j); 
            ff_spin_keys.insert( 
                value_type( s + std::string( "_tran"), rep_pair(*i,*j) )); 
            if( i == j )
              ff_spin_keys.insert( 
                  value_type( s + std::string( "_diag"), rep_pair(*i,*j) )); 
          }

        // snarf some more keys
        ff_allowed_spin_keys 
          = ::radmat::LorentzffFormFactorDecompositionFactoryEnv::all_keys(); 

        // if it is allowed stick it into the map that radmat uses 
        for( it = ff_spin_keys.begin(); it != ff_spin_keys.end(); ++it)
          if( std::find( ff_allowed_spin_keys.begin(), 
                ff_allowed_spin_keys.end(), 
                it->first) != ff_allowed_spin_keys.end() )
            success &= do_reg( it->first, it->second ); 

        registered = true; 
      }

      return success; 
    }

  } // HelicityFormFactorFactoryDecompositionFactoryEnv



} // radmat


#undef CHECK_HELICITY_EXPLICIT
