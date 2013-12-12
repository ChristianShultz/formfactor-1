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
    std::string 
      rotation_group_label(const mom_t &, 
          const mom_t &); 

    typedef Util::SingletonHolder<RotationGroupGenerator<2,4> >
      TheRotationGroupGenerator; 

    bool registerAll();


    rHandle<RotationMatrix_t> 
      get_left_rotation(const mom_t &l, const mom_t &r);

    rHandle<RotationMatrix_t> 
      get_right_rotation(const mom_t &l, const mom_t &r);


  } // LatticeRotationEnv

} // radmat



#endif /* LORENTZFF_CANONICAL_ROTATIONS_H */
