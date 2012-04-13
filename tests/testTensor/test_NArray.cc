// test_NArray.cc -
//
// Monday, March 19 2012
//


#include <complex>
#include "tensor/test_NArray.h"



using namespace tensor;

int 
main(void)
{
  return test_NARRAY<double>() + test_NARRAY<std::complex<double> >();
}



