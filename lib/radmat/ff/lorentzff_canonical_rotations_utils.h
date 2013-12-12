#ifndef LORENTZFF_CANONICAL_ROTATIONS_UTILS_H
#define LORENTZFF_CANONICAL_ROTATIONS_UTILS_H 

#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "itpp/itbase.h"


namespace radmat
{


  // Ax = b
  bool
    check_frame_transformation(const RotationMatrix_t *A,
        const mom_t &x, 
        const mom_t &b);

  // return the transformation from pp to p
  //
  // p = R * pp
  RotationMatrix_t* 
    generate_frame_trasformation(const mom_t &p, const mom_t &pp);


  template<int X, int Y, int Z> 
    mom_t gen_mom(void)
    {
      mom_t p;
      p.resize(3); 
      p[0] = X;
      p[1] = Y; 
      p[2] = Z; 
      return p; 
    }

  int dot_mom_t(const mom_t &left, const mom_t &right);

  mom_t cross_product(const mom_t &l, const mom_t &r);

  bool is_rest(const mom_t &p);

  // test the dot product and volume element are the same
  bool related_by_rotation(const mom_t &left, 
      const mom_t &right, 
      const mom_t &lleft,
      const mom_t &rright,
      bool allow_flip=true);

  itpp::Vec<double> normalize(const mom_t &m); 

}


#endif /* LORENTZFF_CANONICAL_ROTATIONS_UTILS_H */
