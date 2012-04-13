// test_CLONE.cc -
//
// Monday, April  2 2012
//

#include "tensor/tensorbase.h"
#include <iostream>
#include <complex>


using namespace tensor;

int
main(void)
{
  Tensor<double,1> d1,*d2;
  idx_t bd = 4;
  std::vector<idx_t> dim(1,bd);
  d1.create(&dim[0]);

  for(idx_t i = 0; i < bd; ++i)
    d1[i] = i*5;
  
  d2 = d1.clone();

  // std::cout << d1 << std::endl;
  // std::cout << *d2 << std::endl;
  
  delete d2;

  
  Tensor<std::complex<double> , 1> d3,*d4, *d5;
  d3.create(&dim[0]);

  for(idx_t i = 0; i < bd; ++i)
    d3[i] = std::complex<double>(i*5, -i);

  d4 = d3.clone();

  // std::cout << d3 << std::endl;
  // std::cout << *d4 << std::endl;

  d5 = d4->clone();

  delete d4;

  // std::cout << *d5 << std::endl;

  TensorImplBase* base = d5->clone();

  delete d5;

  // std::cout << "check *base = *derived, then cast back and deref" << std::endl; 
  // std::cout << castAndDeref<std::complex<double> , 1>(base) << std::endl;

  delete base;

  std::cout << "All clone tests passed" << std::endl;

  return 0;
}
