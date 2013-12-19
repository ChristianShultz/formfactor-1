/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : test_tens.cc

* Purpose :

* Creation Date : 22-11-2013

* Last Modified : Thu 19 Dec 2013 12:42:52 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/utils/tensor.h"
#include "radmat/utils/handle.h"
#include "radmat/ff/formfactor_factory.h"
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
  
  radmat::FormFactorDecompositionFactoryEnv::registerAll();

  // radmat::RhoRho::RhoRho<-1,-1> rho;
  // std::cout << rho.ff() << std::endl;

  delete pp;

  return 0; 
}

