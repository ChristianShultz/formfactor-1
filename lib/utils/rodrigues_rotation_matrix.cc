// rodrigues_rotation_matrix.cc -
//
// Friday, March 30 2012
//

#include "rodrigues_rotation_matrix.h"
#include "itpp/itbase.h"
#include "pow2assert.h"
#include <math.h>
#include "tensor/tensorbase.h"

using namespace tensor;

// couldn't find it in a real book but wolfram came to the rescue
// http://mathworld.wolfram.com/RodriguesRotationFormula.html

// solves rotated = R*orig for the rotation matrix R where R is a proper rotation
itpp::Mat<double> rodRotMat(const itpp::Vec<double> &orig,
			    const itpp::Vec<double> &rotated)
{
  POW2_ASSERT((orig.size() == rotated.size()) && (orig.size() == 3));
  
  itpp::Vec<double> no = orig/sqrt(orig*orig);
  itpp::Vec<double> nr = rotated/sqrt(rotated*rotated);
  itpp::Vec<double> cross = itpp::cross(no,nr);
  itpp::Mat<double> I3(3,3), cross3(3,3);

  // the general formula will fail for improper rotations, ie parity, we can explicitly make
  // this work by returning -I3 in this case
  bool improper = true;
  double precision = 1e-15;
  for(short i = 0; i < 3; ++i)
    if(fabs(no[i] - nr[i]) > precision)
      improper = false;

  if(improper)
    return -I3;

  double theta = acos(no*nr);

  // need to renormalize this guy to a unit vector
  if(cross*cross > 1e-15)       // numerically compatible with a zero vector -- ~colinear vectors
    cross /= sqrt(cross*cross); // this should never happen with the lattice p vectors

  /*
    std::cout << "no * no = " << no * no << std::endl;
    std::cout << "nr * nr = " << nr * nr << std::endl;
    std::cout << "no * cross = " << no*cross << std::endl;
    std::cout << "nr * cross = " << nr*cross << std::endl;
    std::cout << "cross * cross = " << cross*cross << std::endl;
  */

  I3.zeros();
  cross3.zeros();

  I3(0,0) = 1.;
  I3(1,1) = 1.;
  I3(2,2) = 1.;

  cross3(0,1) = -cross(2);   cross3(1,0) = cross(2);
  cross3(0,2) = cross(1);    cross3(2,0) = -cross(1);
  cross3(1,2) = -cross(0);   cross3(2,1) = cross(0);

  return I3 + cross3*sin(theta) + cross3*cross3*(1.-cos(theta));
}


itpp::Mat<double> rodRotMat(const Tensor<double, 1> &orig, const Tensor<double,1> &rot)
{
  itpp::Vec<double> o(3),r(3);
  for(idx_t i = 0; i < 3; ++i)
    {
      o[i] = orig[i+1];
      r[i] = rot[i+1];
    }
  return rodRotMat(o,r);
}
