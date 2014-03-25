/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : formfactor_kinematic_factors.cc

 * Purpose :

 * Creation Date : 18-03-2014

 * Last Modified : Tue 25 Mar 2014 02:32:47 PM EDT

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
#include <utility>
#include <sstream>
#include <string>
#include <list>

namespace radmat
{


  namespace 
  {

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
        __builtin_trap(); 
        return KF.getRow( tag->gamma_row ); 
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
