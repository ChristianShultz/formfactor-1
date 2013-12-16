#ifndef LORENTZFF_CANONICAL_ROTATIONS_H
#define LORENTZFF_CANONICAL_ROTATIONS_H 


#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "lorentzff_canonical_rotation_groups.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include <complex>
#include <string>
#include "io/adat_xmlio.h"


namespace radmat
{
  namespace LatticeRotationEnv
  {

    // string label of the canonical frame
    std::string 
      rotation_group_label(const mom_t &, 
          const mom_t &); 

    typedef Util::SingletonHolder<RotationGroupGenerator<2,4> >
      TheRotationGroupGenerator; 

    bool registerAll();


    struct FrameOrientation_t
    {
      FrameOrientation_t() {}
      FrameOrientation_t(const mom_t &ccl, const mom_t &ccr, const mom_t &ll, const mom_t &rr)
        : cl(ccl) , cr(ccr) , l(ll) , r(rr)
      {  }
      mom_t cl,cr,l,r; 
    };


    // takes momentum -- returns orientations 
    FrameOrientation_t
      get_frame_orientation(const mom_t &l, const mom_t &r); 




    // these take momentum -- NOT orientations

    // return the rotation back to the reference frame
    //
    //  v_ref = R * v_frame
    //
    rHandle<RotationMatrix_t> 
      get_left_rotation(const mom_t &l, const mom_t &r);

    rHandle<RotationMatrix_t> 
      get_right_rotation(const mom_t &l, const mom_t &r);

    // return the rotation to the canonical direction  
    //
    //  Pcan = R * p_z
    //
    rHandle<RotationMatrix_t>
      get_left_can_frame_rotation(const mom_t &l, const mom_t &r); 
      
    rHandle<RotationMatrix_t>
      get_right_can_frame_rotation(const mom_t &l, const mom_t &r); 

    // for testing only
    //    first pair is can frame
    RotationMatrix_t*
      get_left_rotation(const mom_t &cl, 
        const mom_t &cr, 
        const mom_t &l, 
        const mom_t &r);

    // for testing only
    //    first pair is can frame
    RotationMatrix_t* 
      get_right_rotation(const mom_t &cl, 
        const mom_t &cr, 
        const mom_t &l, 
        const mom_t &r);

  } // LatticeRotationEnv

} // radmat



#endif /* LORENTZFF_CANONICAL_ROTATIONS_H */
