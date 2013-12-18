#ifndef LORENTZFF_CANONICAL_ROTATIONS_UTILS_H
#define LORENTZFF_CANONICAL_ROTATIONS_UTILS_H 

#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "itpp/itbase.h"
#include "hadron/irrep_util.h"


namespace radmat
{


  // Ax = b
  bool
    check_frame_transformation(const RotationMatrix_t *A,
        const mom_t &x, 
        const mom_t &b,
        bool print_on_false=false);

  //
  //  (R ll ?= l) && (R rr ?= r)
  //
  bool check_total_frame_transformation(const RotationMatrix_t *R, 
      const mom_t &l,
      const mom_t &r,
      const mom_t &ll, 
      const mom_t &rr,
      bool print_on_false=false); 


  // are they colinear
  bool colinear_momentum(const mom_t &l, const mom_t &r); 

  // return the transformation from pp to p
  //
  // first = R * second
  //
  // R is generated from the lg rotations in adat as 
  // R = R_lat(first) * R^-1_lat(second)
  RotationMatrix_t* 
    generate_frame_transformation(const mom_t &first, const mom_t &second);

  //  TRIaxial Attitude Determination 
  //
  // fails if @ rest 
  //
  // returns the matrix that solves 
  // 
  //     l = R * ll 
  //     r = R * rr
  //
  // NB: if the momentum are colinear then 
  //     they are in the same lg so 
  //     return
  //      generate_frame_transformation(l,ll)
  RotationMatrix_t* 
    generate_triad_rotation_matrix(const mom_t &l, const mom_t &r, 
        const mom_t &ll, const mom_t &rr); 

  // parameterize the triad rotaion in terms
  // of z-y-z euler angles 
  Hadron::CubicCanonicalRotation_t
    generate_euler_angles(const RotationMatrix_t*); 

  // kill things below a thresh 
  void clean_up_rot_mat(RotationMatrix_t *);

  // useful 
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

  // integer dot product
  int dot_mom_t(const mom_t &left, const mom_t &right);

  // integer cross product
  mom_t cross_product(const mom_t &l, const mom_t &r);

  // self explanatory 
  bool is_rest(const mom_t &p);

  // throw if passed a rest momentum 
  bool related_by_rotation(const mom_t &left, 
      const mom_t &right, 
      const mom_t &lleft,
      const mom_t &rright,
      bool allow_flip=true);

  // return a unit vector in the direction of m
  itpp::Vec<double> normalize(const mom_t &m); 

  // take a determinant to check it it is +1
  double determinant(const RotationMatrix_t* ); 

}


#endif /* LORENTZFF_CANONICAL_ROTATIONS_UTILS_H */
