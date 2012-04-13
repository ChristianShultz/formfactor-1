// test_Tensor.cc -
//
// Monday, March 26 2012
//


#include "tensor/test_Tensor.h"

#include <iostream>
#include <complex>

using namespace tensor;


int 
main(void)
{
  return test_tensor<double>() + test_tensor<std::complex<double> >();
}
