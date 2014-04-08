/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_kinematic_factors.cc

 * Purpose :

 * Creation Date : 18-03-2014

 * Last Modified : Tue 08 Apr 2014 10:14:13 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "formfactor_kinematic_factors.h"
#include "formfactor_factory.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble/semble_meta.h"
#include "radmat/utils/printer.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/redstar_interface/redstar_cartesian_interface.h"
#include "radmat/redstar_interface/redstar_photon_props.h"
#include <utility>
#include <sstream>
#include <string>
#include <list>

namespace radmat
{


  namespace 
  {

    struct three_point_tag_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "three_point_tag_printer " << msg << std::endl;}
    };

    struct subduce_info_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "subduce_info_printer " << msg << std::endl;}
    };

    struct kgen_info_printer
    {
      static void print(const std::string &msg)
      {}
      // { std::cout << "kgen_info_printer " << msg << std::endl;}
    };

    struct helicity_printer
    {
      static void print(const std::string &msg) 
      {}
      // { std::cout << "helicity_printer " << msg << std::endl; }
    }; 

    std::string to_string(const int i)
    {std::stringstream ss; ss << i; return ss.str();}

    std::string to_string(const std::complex<double> &cd)
    {std::stringstream ss; ss << cd; return ss.str();}

    std::string to_string(const itpp::Mat<std::complex<double> > &m )
    {std::stringstream ss; ss << m; return ss.str(); }


    FFKinematicFactors_t::KinematicFactorRow
      handle_subduce_J0( const FFKinematicFactors_t::KinematicFactorMatrix &KF , 
          const int row, 
          const ADATXML::Array<int> &q, 
          const std::string &sub_key)
      {
        const SubduceTableMap::irrep_sub_table * table; 
        table = TheSmarterSubduceTableMap::Instance().get_table(sub_key); 
        SubduceTableMap::sub_list subduce = table->query_1( row ); 

        // subduce is a list of pairs of stl complex and int  

        // this must be a 1 elem list
        POW2_ASSERT(subduce.size() == 1);

        // the scalar part lives in the 0 slot
        return KF.getRow(0);  
      }


    // rows are + 0 - 
    FFKinematicFactors_t::KinematicFactorMatrix 
      move_to_helicity( const FFKinematicFactors_t::KinematicFactorMatrix &KF , 
          const ADATXML::Array<int> &q)
      {
        FFKinematicFactors_t::KinematicFactorMatrix ret; 

        // single precision zero  
        double tolerance = 1e-6; 

        // rows are + 0 -  , cols are x y z 
        itpp::Mat<std::complex<double> > eps = itpp::round_to_zero(eps3d(q,PHOTON_CREATE),tolerance); 


        // intermediate variables 
        FFKinematicFactors_t::KinematicFactorRow p,z,m; 
        p = KF.getRow(0); 
        p.zeros(); 
        z = p;
        m = p; 

        // ensems are big, avoid a copy if the coeff is zero-ish 
        std::complex<double> complex_zero(0.,0.); 

        // take the linear combinations, 
        // bung them into named vectors
        //
        // NB: KF is a 4 tensor, eps is a 3 tensor
        for(int i = 0; i < 3; ++i)
        {
          if( eps(0,i) != complex_zero )
            p += eps(0,i) * KF.getRow(i+1); 
          if( eps(1,i) != complex_zero )
            z += eps(1,i) * KF.getRow(i+1); 
          if( eps(2,i) != complex_zero )
            m += eps(2,i) * KF.getRow(i+1); 
        }


        // dump the vectors into the return matrix 
        ret.reDim(p.getB(), 3, p.getN()); 
        for(int j = 0; j < p.getN(); ++j)
        {
          ret.loadEnsemElement(0,j,p.getEnsemElement(j));
          ret.loadEnsemElement(1,j,z.getEnsemElement(j));
          ret.loadEnsemElement(2,j,m.getEnsemElement(j));
        }

        //  printer_function<helicity_printer>( "\n" +  to_string(eps) ) ; 
        //  printer_function<helicity_printer>( " - KF \n" + to_string( KF.mean() ) ); 
        //  printer_function<helicity_printer>( " - eps * KF\n" + to_string( ret.mean() ) ); 


        return ret; 
      }

    FFKinematicFactors_t::KinematicFactorRow
      handle_subduce_J1( const FFKinematicFactors_t::KinematicFactorMatrix &KF , 
          const int row, 
          const ADATXML::Array<int> &q, 
          const std::string &sub_key)
      {
        FFKinematicFactors_t::KinematicFactorRow ret; 

        // subduce is a list of pairs of stl complex and int  
        const SubduceTableMap::irrep_sub_table * table; 
        table = TheSmarterSubduceTableMap::Instance().get_table(sub_key); 
        SubduceTableMap::sub_list subduce = table->query_1( row ); 

        // + -> row 0, 0 -> row 1, - -> row 2 of the matric Hel
        FFKinematicFactors_t::KinematicFactorMatrix Hel = move_to_helicity( KF, q); 

        // init it 
        ret = Hel.getRow(0); 
        ret.zeros(); 

        std::complex<double> complex_zero(0.,0.); 

        // loop the coeffs and put it together 
        SubduceTableMap::sub_list::const_iterator it; 
        for(it = subduce.begin(); it != subduce.end(); ++it)
          if( it->first != complex_zero )
          {
            printer_function<subduce_info_printer>( 
                to_string(it->first) + " x [" + to_string(it->second) +"]");  
            ret += it->first * Hel.getRow(1 - it->second); 
          }

        return ret; 
      }


    FFKinematicFactors_t::KinematicFactorRow
      handle_three_point_lorentz_insertion( 
          const FFKinematicFactors_t::KinematicFactorMatrix &KF,
          const ThreePointDataTag * tag)
      {
        return KF.getRow( tag->gamma_row ); 
      }

    FFKinematicFactors_t::KinematicFactorRow
      handle_three_point_cubic_insertion( 
          const FFKinematicFactors_t::KinematicFactorMatrix &KF,
          const ThreePointDataTag * tag, 
          const rHandle<Rep_p> & gamma)
      {
        // return variable
        FFKinematicFactors_t::KinematicFactorRow ret = KF.getRow(0); 
        ret.zeros(); 

        printer_function<three_point_tag_printer>(gamma->rep_type()); 

        // create the key  -- cubic part
        const CubicRep * g_rep_ptr; 
        g_rep_ptr = dynamic_cast<const CubicRep*>(gamma.get_ptr()); 
        rHandle<CubicRep> g_rep( g_rep_ptr->clone() ); 

        // create the key  -- lorentz part
        rHandle<Rep_p> g_origin_prim = tag->origin_rep.gamma(); 
        const LorentzRep * g_origin_rep; 
        g_origin_rep = dynamic_cast<const LorentzRep*>(g_origin_prim.get_ptr()); 
        rHandle<LorentzRep> g_origin( g_origin_rep->clone() ); 

        // subduction table key 
        std::string key = make_subduce_table_map_id(g_origin,g_rep);  

        printer_function<three_point_tag_printer>( "subduce key->" + key ); 
        int g_row = tag->gamma_row; 
        ADATXML::Array<int> q = tag->q; 

        if( g_origin->rep_spin() == 0 )
        {
          ret = handle_subduce_J0( KF, g_row, q, key); 
        }
        else if( g_origin->rep_spin() == 1 )
        {
          ret = handle_subduce_J1( KF, g_row, q, key); 
        }
        else
        {
          // die with a seg fault, photons are spin 0 or 1 for me,
          // I dont care about you
          std::cout << __PRETTY_FUNCTION__ << " err: " 
            << g_origin->rep_id() << " is not a photon" << std::endl;
          __builtin_trap(); 
        }
        return ret; 
      }

    FFKinematicFactors_t::KinematicFactorRow
      handle_three_point_insertion( const FFKinematicFactors_t::KinematicFactorMatrix &KF,
          const ThreePointDataTag * tag)
      {
        rHandle<Rep_p> gamma = tag->data_rep.gamma(); 
        if( gamma->rep_type() == Stringify<CubicRep_t>() )
          return handle_three_point_cubic_insertion( KF, tag, gamma); 
        else
          return handle_three_point_lorentz_insertion( KF , tag); 
      }


    FFKinematicFactors_t::KinematicFactorRow 
      generate_three_point_factors( const DataTagPrimitive *ptr)
      {
        const ThreePointDataTag * tag = dynamic_cast<const ThreePointDataTag*>(ptr); 
        rHandle<FormFactorBase_t> mat_elem; 
        mat_elem = FormFactorDecompositionFactoryEnv::callFactory( tag->mat_elem_id); 

        // size params 
        int nfacs = mat_elem->nFacs(); 
        int nbins = tag->left_E.size(); 
        FFKinematicFactors_t::KinematicFactorMatrix KF(nbins,4,nfacs); 
        KF.zeros(); 

        // scale down the energies, momentum have zero variance 
        //   NB: need to use assignment here
        ENSEM::EnsemReal left_E , right_E; 
        left_E = ENSEM::rescaleEnsemDown( tag->left_E ); 
        right_E = ENSEM::rescaleEnsemDown( tag->right_E ); 

        // mom and row params 
        ADATXML::Array<double> left_mom(3),right_mom(3); 
        int left_row, right_row; 
        left_row = tag->left_row; 
        right_row = tag->right_row; 

        // the momentum unit
        double mom_factor = tag->mom_fac; 

        for(int i =0; i < 3; ++i )
        {
          left_mom[i] = mom_factor * tag->left_mom[i]; 
          right_mom[i] = mom_factor * tag->right_mom[i]; 
        }

        // this also checks size of righty vs size of lefty
        POW2_ASSERT_DEBUG( (nbins == right_E.size()) && (nfacs > 0) );

        printer_function<kgen_info_printer>(tag->mat_elem_id); 
        printer_function<kgen_info_printer>("left_row " + to_string(left_row)); 
        printer_function<kgen_info_printer>("right_row " + to_string(right_row)); 

        // loop over cfgs and use the wrapper to operator()
        for(int bin = 0; bin < nbins; bin++)
          KF[bin] = mat_elem->wrapper( left_E.elem(bin) , left_mom , left_row,
              right_E.elem(bin) , right_mom , right_row , mom_factor);

        // scale up return data
        KF.rescaleSembleUp();

        // now handle the representation of the photon 
        //     remember that all decomps return a lorentz
        //     vector 

        return handle_three_point_insertion( KF, tag ); 
      }





  }


  // the tag contains all knowledge ( well at least of the radmat related variety ) 
  FFKinematicFactors_t::KinematicFactorRow 
    FFKinematicFactors_t::genFactors(const DataTagPrimitive * ptr)
    {

      typedef KinematicFactorRow (*ptr_type)(const DataTagPrimitive *); 
      std::map<std::string, ptr_type> options_map; 

      // may someday try nucleon types or some fancy luscher garbage?
      options_map.insert( std::make_pair( Stringify<ThreePointDataTag>(), 
            &generate_three_point_factors ) ) ; 

      std::map<std::string,ptr_type>::const_iterator it; 
      it = options_map.find( ptr->type() ); 

      if( it == options_map.end() )
      {
        std::cout << __PRETTY_FUNCTION__ << ": error," 
          << " I dont know what a " << ptr->type() 
          << " is, options were " << std::endl;

        for(it = options_map.begin(); it != options_map.end(); ++it)
          std::cout << it->first << std::endl; 

        exit(1); 
      }

      ptr_type foo = it->second; 

      return (*foo)(ptr); 
    }



} // radmat 
