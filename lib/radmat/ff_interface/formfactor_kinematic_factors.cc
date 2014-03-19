/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : formfactor_kinematic_factors.cc

* Purpose :

* Creation Date : 18-03-2014

* Last Modified : Tue 18 Mar 2014 03:00:26 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "formfactor_kinematic_factors.h"
#include "photon_subduction.h"

namespace radmat
{

  FFKinematicFactors_t::KinematicFactorRow
    FFKinematicFactors_t::gamma_rep( const FFKinematicFactors_t::KinematicFactorMatrix &M,
        const std::string &rep, 
        const std::string &spher_rep,
        const int row,
        const ADATXML::Array<int> &q) const
    {
      FFKinematicFactors_t::KinematicFactorRow ret; 
      ret = M.getRow(0); 
      ret.zeros(); 

      lorentz_to_cubic_t sub = photon_subduction(rep,spher_rep,row,q);

      lorentz_to_cubic_t::const_iterator it; 
      for(it = sub.begin(); it != sub.end(); ++it)
        ret += it->m_coeff * M.getRow(it->m_obj); 

      return ret; 
    }


} // radmat 
