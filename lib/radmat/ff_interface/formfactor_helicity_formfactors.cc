/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_helicity_formfactors.cc

 * Purpose :

 * Creation Date : 22-02-2014

 * Last Modified : Sat 22 Feb 2014 04:53:09 PM EST

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
#include <sstream>
#include <algorithm>

#define CHECK_HELICITY_EXPLICIT

namespace radmat
{

  // sanity
#ifdef CHECK_HELICITY_EXPLICIT
  namespace
  {
    bool check_helicity(const int h, const int row, const int J)
    {
      return h == J - row + 1; 
    }
  }
#endif

  // use whatever is in the lorentzff_canonical_cont_spin decomp
  itpp::Mat<std::complex<double> >
    HelicityFormFactor::generate_ffs( const MomRowPair_t &lefty, 
        const MomRowPair_t &righty, 
        const double mom_fac) const
    {
      FormFactorBase_t::recipe_h recipe = FormFactorBase_t::get_recipe(); 
      POW2_ASSERT( recipe->id() == Stringify<HelicityFormFactorRecipe_t>() ); 
      const HelicityFormFactorRecipe_t * r;
      r = dynamic_cast< const HelicityFormFactorRecipe_t* >( recipe.get_ptr() );

#ifdef CHECK_HELICITY_EXPLICIT
      check_helicity( lefty.second, r->left_row() , r->left_spin() ); 
      check_helicity( righty.second, r->right_row(), r->right_spin() ); 
#endif

      return r->mat->operator()(lefty,righty,mom_fac); 
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
      bool do_reg( const std::string &s, const rep_pair &r)
      {
        return true; 
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
          {
            std::cout << "found " << it->first << std::endl;
            success &= do_reg( it->first, it->second ); 
          }

        registered = true; 
      }

      return success; 
    }

  } // HelicityFormFactorFactoryDecompositionFactoryEnv



} // radmat


#undef CHECK_HELICITY_EXPLICIT
