/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : rotation_group_generator.cc

* Purpose :

* Creation Date : 14-04-2014

* Last Modified : Wed 16 Apr 2014 09:17:43 AM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "rotation_group_generator.h"
#include "radmat/utils/utils.h"



namespace radmat
{

  std::ostream & operator<<(std::ostream &o, const mom_key &k)
  { return o << k.x << k.y << k.z ;}

  std::ostream & operator<<(std::ostream &o, const mom_pair_key &k)
  { return o << "l" << k.l << " r" << k.r; }


  namespace LatticeRotationEnv
  {

    namespace
    {
      bool local_registration = false; 
    }


    bool registerAll(void)
    {
      printer_function<debug_print_reg_all>(__PRETTY_FUNCTION__); 
      bool success = true; 

      if ( !!! local_registration)
      {
        success &= TheRotationGroupGenerator::Instance().registerAll(); 
        local_registration = true; 
      }

      if( !!! success )
      {
        throw std::string("reg error in LatticeRotationEnv");
      }

      return success; 
    }

    //////////////////////////////////////////////////////////////////////

    std::pair<mom_t,mom_t>
      rotation_group_key(const mom_t &l, const mom_t &r)
      {
        return TheRotationGroupGenerator::Instance().canonical_frame(l,r);  
      }
  
    std::string 
      rotation_group_label(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> k = rotation_group_key(l,r); 
        std::stringstream ss; 
        ss << "lefty" << k.first[0] << k.first[1] << k.first[2] 
          << "righty" << k.second[0] << k.second[1] << k.second[2]; 
        return ss.str(); 
      }

  } // LatticeRotationEnv

} // radmat 

