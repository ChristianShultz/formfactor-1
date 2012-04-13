// test_rodrigues.cc -
//
// Friday, March 30 2012
//

#include <iostream>
#include "itpp/itbase.h"
#include "utils/rodrigues_rotation_matrix.h"
#include "utils/pow2assert.h"


int 
main(void)
{
  itpp::Vec<double> orig(3), rot(3);

  // test a random rotation
  orig = itpp::randu(3);
  orig = orig/sqrt(orig*orig);

  rot = itpp::randu(3);
  rot = rot/sqrt(rot*rot);

  itpp::Mat<double> R = rodRotMat(orig,rot);

  for(int i = 0; i < 3; ++i)
    POW2_ASSERT(fabs(rot(i) - (R*orig)(i)) < 1e-15);


  // test a rotation by pi
  rot = -orig;
  R = rodRotMat(orig,rot);
  
  std::cout << "Rodrigues Rotation Matrix test successful" << std::endl;

  return 0;
}
