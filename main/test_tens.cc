/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : test_tens.cc

* Purpose :

* Creation Date : 22-11-2013

* Last Modified : Fri 22 Nov 2013 02:03:48 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include <iostream>

int 
main(void)
{
 
  radmat::Tensor<double,1> foo((radmat::TensorShape<1>())[4],0.);  
  radmat::Tensor<double,1> bar((radmat::TensorShape<1>())[4],1.);  
  radmat::Tensor<double,1> baz;
  
  baz = foo + bar; 

  radmat::rHandle<radmat::Tensor<double,1> > p ( foo.clone() ); 
  radmat::Tensor<double,1> *pp = baz.clone(); 

  delete pp;

  return 0; 
}

