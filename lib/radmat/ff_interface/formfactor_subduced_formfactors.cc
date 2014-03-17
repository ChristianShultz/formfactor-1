/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_subduced_formfactors.cc

 * Purpose :

 * Creation Date : 13-03-2014

 * Last Modified : Mon 17 Mar 2014 01:19:41 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "formfactor_subduced_formfactors.h"
#include "formfactor.h"
#include "formfactor_helicity_formfactors.h"
#include "formfactor_factory.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/printer.h"
#include "hadron/subduce_tables_factory.h"
#include "semble/semble_meta.h"
#include <sstream>
#include <algorithm>


typedef radmat::TheFormFactorInjectionRecipeCookbook Cookbook; 
typedef radmat::TheFormFactorInjectionRecipeFactory Factory; 



namespace radmat
{

  namespace 
  {
    
    struct subduce_table_printer
    {
      static void print(const std::string &msg)
      { std::cout << "subduce_table_printer " + msg << std::endl; }
    };


    const SubduceTableMap::irrep_sub_table* 
      query_subduce_table(const std::string &table_id)
      {
        SubduceTableMap::map_t::const_iterator it; 
        it = TheSmarterSubduceTableMap::Instance().mappy.find(table_id);
        if( it == TheSmarterSubduceTableMap::Instance().mappy.end() )
        {
          printer_function<console_print>( "missed subduce table " + table_id ); 
          throw std::string("missed subduce table"); 
        }

        return it->second; 
      }

    SubduceTableMap::irrep_sub_table*
    subduce_rest(const std::string &table_id, 
        const rHandle<SpherRep_p> &cont, 
        const rHandle<CubicRep_p> &cub)
    {
      SubduceTableMap::irrep_sub_table::sub_table tab; 
      Hadron::SubduceTable* sub = Hadron::TheSubduceTableFactory::Instance().createObject(table_id); 
      std::complex<double> zero(0.,0.);

      tab.resize( cub->dim() ); 

      for(int row = 0; row < cub->dim(); ++row )
      {
        SubduceTableMap::sub_list entry; 
        for(int lambda = 0; lambda < 2*cont->rep_spin() + 1; ++lambda)
        {
          // subduce tables are FORTRAN based 
          std::complex<double> weight = SEMBLE::toScalar( (*sub)(row+1,lambda+1) ); 
          if(weight == zero)
            continue; 

          entry.push_back(std::make_pair(weight,cont->rep_spin() - lambda) ); 
        }
        tab[row] = entry; 
      }

      delete sub; 

      SubduceTableMap::irrep_sub_table* ret = new SubduceTableMap::irrep_sub_table(tab,cont,cub); 
      return ret; 
    }

    SubduceTableMap::irrep_sub_table* 
      subduce_flight(const std::string &table_id, 
        const rHandle<SpherRep_p> &cont, 
        const rHandle<CubicRep_p> &cub)
      { 
      SubduceTableMap::irrep_sub_table::sub_table tab; 
      Hadron::SubduceTable* sub = Hadron::TheSubduceTableFactory::Instance().createObject(table_id); 
      const CubicRepFlight_p * cubb;
      cubb = dynamic_cast< const CubicRepFlight_p *> ( cub.get_ptr() ); 
      std::complex<double> zero(0.,0.);
      double eta_p = double(cont->rep_eta_p()); 
      tab.resize( cub->dim() ); 
      int helicity = cubb->helicity();


      if ( helicity == 0 )
      {
        SubduceTableMap::sub_list entry;
        std::complex<double> weight = SEMBLE::toScalar( (*sub)(1,1) );
        entry.push_back( std::make_pair( weight, helicity)); 
        tab[0] = entry;
      }
      else
      {
        for(int row = 0; row < cub->dim(); ++row)
        {
          SubduceTableMap::sub_list entry;
          std::complex<double> weightp = SEMBLE::toScalar( (*sub)(row+1,1));
          std::complex<double> weightm = SEMBLE::toScalar( (*sub)(row+1,2));
          
          entry.push_back(std::make_pair(weightp,helicity));
          entry.push_back(std::make_pair( eta_p * weightm , -helicity ) ); 

          tab[row] = entry; 
        }
      }

      delete sub;

      SubduceTableMap::irrep_sub_table* ret = new SubduceTableMap::irrep_sub_table(tab,cont,cub); 
      return ret; 
    }

    std::string make_subduce_table_map_id(const rHandle<SpherRep_p> &cont, 
        const rHandle<CubicRep_p> &cub)
    {
      return cont->rep_id() + cub->rep_id(); 
    }

    void add_subduce_table( const std::string &table_id, 
        const rHandle<SpherRep_p> &cont, 
        const rHandle<CubicRep_p> &cub)
    {
      std::string map_id = make_subduce_table_map_id(cont,cub); 

      if( TheSmarterSubduceTableMap::Instance().mappy.find(map_id)
          == TheSmarterSubduceTableMap::Instance().mappy.end() )
      {
        printer_function<subduce_table_printer>( "made a " + map_id ); 

        if( cub->rep_g() == Stringify<Oh>() )
          TheSmarterSubduceTableMap::Instance().mappy.insert(
              std::make_pair( map_id , subduce_rest(table_id, cont , cub ) ) ); 
        else
          TheSmarterSubduceTableMap::Instance().mappy.insert(
              std::make_pair( map_id , subduce_flight(table_id, cont , cub ) ) ); 
      }
    }

  } // anonomyous 

  // all the nasty comes together here 
  itpp::Mat<std::complex<double> >
    SubducedFormFactor::generate_ffs(const MomRowPair_t &lefty, 
        const MomRowPair_t &righty, 
        const double mom_fac) const 
    {
      FormFactorBase_t::recipe_h recipe_h = FormFactorBase_t::get_recipe(); 
      POW2_ASSERT( recipe_h->id() == Stringify<SubducedFormFactorRecipe_t>() ); 
      const SubducedFormFactorRecipe_t *recipe; 
      recipe = dynamic_cast<const SubducedFormFactorRecipe_t *>(recipe_h.get_ptr()); 

      // this is the return object, it is a 4 x nff matrix, zero it
      // out explicitly and then just sum on the heicities with 
      // a weight to generate a subduced decomposition 
      itpp::Mat<std::complex<double> > ret(4,recipe->nFacs());
      ret.zeros(); 

      // the rows that we store internally are 0 based 
      int left_row = lefty.second - 1; 
      int right_row = righty.second - 1; 
      std::complex<double> cmplx_zero(0.,0.); 

      // left and right representations to make the strings to 
      // call the factory 
      rHandle<CubicRep_p> l_cub_rep , r_cub_rep;
      rHandle<SpherRep_p> l_sph_rep , r_sph_rep; 

      l_cub_rep = recipe->lefty;
      r_cub_rep = recipe->righty;
      l_sph_rep = recipe->hel.lefty;
      r_sph_rep = recipe->hel.righty; 

      // generate the id for the subduce table map 
      std::string l_map_id , r_map_id; 
      l_map_id = make_subduce_table_map_id( l_sph_rep , l_cub_rep );
      r_map_id = make_subduce_table_map_id( r_sph_rep , r_cub_rep ); 

      // the subduction tables -- these were filled out 
      // earlier in the registration stage 
      const SubduceTableMap::irrep_sub_table* l_table;
      const SubduceTableMap::irrep_sub_table* r_table;

      // find the table we want 
      l_table = query_subduce_table(l_map_id);
      r_table = query_subduce_table(r_map_id); 

      // got a const ptr, use a const iterator 
      SubduceTableMap::sub_list::const_iterator l,r; 

      // make the linear combinations 
      //     -- this call eventually descends all the way down to the 
      //     lorentzff_canonical classes, these classes know about how
      //     to transform under rotations so all the hard work is hidden 
      for( l = l_table->sub[left_row].begin(); l != l_table->sub[left_row].end(); ++l)
        for( r = r_table->sub[right_row].begin(); r != r_table->sub[right_row].end(); ++r)
         ret += ( ( std::conj(l->first) * r->first ) 
                * recipe->hel.mat->operator()( std::make_pair( lefty.first , l->second),
                      std::make_pair( righty.first , r->second ),
                      mom_fac)); 

      return ret; 
    }




  namespace SubducedFormFactorDecompositionFactoryEnv
  {

    namespace
    {
      struct helicity_recipe_printer
      {
        static void print(const std::string &s)
        {}
      //  { std::cout << "subduced form factors, found a " << s << std::endl;}
      };

      struct subduce_printer
      {
        static void print(const std::string &s)
        {}
      //  { std::cout << "found subduction table " << s << std::endl;}
      };

      struct subduce_reg_printer
      {
        static void print(const std::string &s)
        {}
      //  { std::cout << "subduced form factors, regging  " << s << std::endl;}
      };


      // loop the cookbook and pull out helicity recipes for subduction 
      std::vector<rHandle<HelicityFormFactorRecipe_t> > 
        get_all_helicity_recipies()
        {
          std::vector<rHandle<HelicityFormFactorRecipe_t> > ret;
          PtrRecipeHolder::map_t::iterator r; 
          for( r = Cookbook::Instance().mappy.begin();
              r != Cookbook::Instance().mappy.end();
              ++r)
          {
            if( r->second->id() == ::radmat::Stringify<HelicityFormFactorRecipe_t>() )
            {
              printer_function<helicity_recipe_printer>(r->first); 
              HelicityFormFactorRecipe_t *ptr = dynamic_cast<HelicityFormFactorRecipe_t*>(r->second);
              ret.push_back( 
                  rHandle<HelicityFormFactorRecipe_t>( new HelicityFormFactorRecipe_t( *ptr ) ) ); 
            }
          }

          return ret; 
        }

      std::vector<std::string> 
        all_subduce_table_keys()
        {
          return Hadron::TheSubduceTableFactory::Instance().keys(); 
        }

      std::vector<std::string> 
        spher_strings(const rHandle<SpherRep_p> &s)
        {
          std::vector<std::string> ret; 
          std::stringstream ss; 
          ss << "J" << s->rep_spin(); 
          ret.push_back(ss.str());
          ss.str("");

          ss << "H0"; 
          s->rep_eta_p() == 1 ? ss << "+" : ss << "-";
          ret.push_back(ss.str());


          for(int i = 1; i <= s->rep_spin(); ++i)
          {
            ss.str("");
            ss << "H" << i; 
            ret.push_back(ss.str()); 
          }

          return ret; 
        };

      // string is the subduce table key, rep is what radmat calls the rep 
      //
      // do a double loop over possible left and right part of subduce table
      // keys to see if it exists 
      std::vector<std::pair<std::string,rHandle<CubicRep_p> > >
        possible_subductions( const rHandle<SpherRep_p> &s )
        {
          std::vector<std::pair<std::string,rHandle<CubicRep_p> > > ret; 
          std::string delim = "->";
          std::string tail = ",1";
          std::vector<std::string> subduce_tables = all_subduce_table_keys();

          std::vector<std::string> sph_strs = spher_strings(s);
          std::vector<std::string> cub_strs = ::radmat::CubicInvariantsFactoryEnv::all_keys(); 
          std::vector<std::string>::const_iterator sph_it,cub_it; 

          for(sph_it = sph_strs.begin(); sph_it != sph_strs.end(); ++sph_it)
            for(cub_it = cub_strs.begin(); cub_it != cub_strs.end(); ++cub_it)
            {
              std::string id = *sph_it + delim + *cub_it + tail; 
              if( std::find( subduce_tables.begin(), subduce_tables.end(), id ) != subduce_tables.end() )
              {
                printer_function<subduce_printer>(id); 
                ret.push_back(std::make_pair(id,::radmat::CubicInvariantsFactoryEnv::callFactory(*cub_it)));
              }
            }

          return ret; 
        }

      FormFactorRecipe_t* gen_recipe( const HelicityFormFactorRecipe_t &h_rep, 
          const rHandle<CubicRep_p> &lefty,
          const rHandle<CubicRep_p> &righty,
          const std::string &left_table, 
          const std::string &right_table)
      {
        add_subduce_table( left_table , h_rep.lefty , lefty );
        add_subduce_table( right_table , h_rep.righty, righty); 
        return new SubducedFormFactorRecipe_t(h_rep,lefty,righty,left_table,right_table); 
      }

      FFAbsBase_t * callback(const std::string &recipe_id)
      {
        PtrRecipeHolder::map_t::iterator r; 
        r = Cookbook::Instance().mappy.find(recipe_id); 

        if( r == Cookbook::Instance().mappy.end())
        {
          printer_function<console_print>( "missed " + recipe_id ); 
          throw std::string("recipe error");
        }

        if( r->second->id() != Stringify<SubducedFormFactorRecipe_t>())
        {
          printer_function<console_print>( "expected a " + Stringify<SubducedFormFactorRecipe_t>());
          printer_function<console_print>( "got a " + r->second->id() ); 
          throw std::string("recipe error");  
        }

        SubducedFormFactorRecipe_t * ptr = dynamic_cast<SubducedFormFactorRecipe_t*>(r->second);
        rHandle<SubducedFormFactorRecipe_t> recipe( new SubducedFormFactorRecipe_t( *ptr ) ); 

        return new SubducedFormFactor( recipe ); 
      }



      bool do_reg( const rHandle<HelicityFormFactorRecipe_t> &h_rep, 
          const std::pair<std::string,rHandle<CubicRep_p> > &lefty,
          const std::pair<std::string,rHandle<CubicRep_p> > &righty)
      {
        std::string reg_id = h_rep->reg_id() + "__" + build_id(lefty.second,righty.second); 
        printer_function<subduce_reg_printer>(reg_id); 
        Cookbook::Instance().mappy.insert(
            std::make_pair(reg_id,
              gen_recipe( *(h_rep.get_ptr()),
                lefty.second,
                righty.second,
                lefty.first,
                righty.first))); 

        bool b = Factory::Instance().registerObject(reg_id,callback);
        if( !!! b )
          printer_function<console_print>( reg_id + " failed to register" ); 

        return b; 
      }

    } // anonomyous 


    std::string build_id(const rHandle<CubicRep_p> &lefty, const rHandle<CubicRep_p> &righty)
    {
      return lefty->rep_id() + "," + righty->rep_id(); 
    }

    namespace
    {
      bool local_registration = false; 
    } 


    // pump the factory 
    bool registerAll()
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if( !!! local_registration )
      {
        std::vector<rHandle<HelicityFormFactorRecipe_t> > h_recipies = get_all_helicity_recipies();
        std::vector<rHandle<HelicityFormFactorRecipe_t> >::const_iterator h_recipe_it; 

        for( h_recipe_it = h_recipies.begin(); h_recipe_it != h_recipies.end(); ++h_recipe_it)
        {
          std::vector<std::pair<std::string,rHandle<CubicRep_p> > > lefty, righty;
          std::vector<std::pair<std::string,rHandle<CubicRep_p> > >::const_iterator l, r;

          lefty = possible_subductions( (*h_recipe_it)->lefty );
          righty = possible_subductions( (*h_recipe_it)->righty );

          for( l = lefty.begin(); l != lefty.end(); ++l)
            for(r = righty.begin(); r != righty.end(); ++r)
              success &= do_reg(*h_recipe_it , *l, *r ); 

        } // h_recipe_it 
      } // local_registration 

      return success; 
    }


  } // SubducedFormFactorDecompositionFactoryEnv

} // radmat

