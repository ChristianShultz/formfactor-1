// test_Print.cc -
//
// Friday, March 30 2012
//

#include "radmat/tensor/tensorbase.h"

using namespace radmat;

int
main(void)
{
  Tensor<double,1> dum;

  std::cout << dum << std::endl;

  std::vector<idx_t> dim(1,4);

  dum.create(&dim[0]);


  std::cout << dum << std::endl;

  return 0;
}
