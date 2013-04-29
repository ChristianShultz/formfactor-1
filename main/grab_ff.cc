/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : grab_ff.cc

* Purpose :

* Creation Date : 22-04-2013

* Last Modified : Mon Apr 22 14:04:01 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>

#include "radmat/ff/formfactor_factory.h"
#include "radmat/utils/tensor.h"


using namespace radmat;
using namespace std; 

struct inp
{
  string matElemID; 
  Tensor<double,1> pf,pi;
  double pfac;
};


inp usage(int argc, char *argv[])
{
  if(argc != 11) 
  {
    std::cerr << "usage: grab ff <matElemID> <Ef0> <nx ny nz>  <Ei0> <nx ny, nz> <2pi/xi/nspace>" << std::endl;
    exit(1); 
  }


  inp foo;

  { std::stringstream ss(argv[1]); ss >> foo.matElemID; }
  { std::stringstream z(argv[2]),o(argv[3]),t(argv[4]),th(argv[5]);
    Tensor<double,1> fred((TensorShape<1>())[4],0.);
    z >> fred[0];
    o >> fred[1];
    t >> fred[2];
    th >> fred[3];
     
    foo.pf = fred; }

    { std::stringstream z(argv[6]),o(argv[7]),t(argv[8]),th(argv[9]);
    Tensor<double,1> fred((TensorShape<1>())[4],0.);
    z >> fred[0];
    o >> fred[1];
    t >> fred[2];
    th >> fred[3];
     
    foo.pi = fred; }

    {std::stringstream ss(argv[10]); ss >> foo.pfac;}


    for(int i = 1; i < 4; ++i)
    {
      foo.pf[i] *= foo.pfac;
      foo.pi[i] *= foo.pfac;
    }

    return foo; 
}


int main(int argc, char *argv[])
{

  inp fred = usage(argc, argv); 
  
  FormFactorDecompositionFactoryEnv::registerAll(); 

  ADAT::Handle<ffBase_t<std::complex<double> > > foo = FormFactorDecompositionFactoryEnv::callFactory(fred.matElemID); 

  std::cout << " The list of ffs is : \n" << foo->ff() << std::endl;

  std::cout << "by value (row index is lorentz, col is FF num) \n" << (*foo)(fred.pf,fred.pi,fred.pfac);

  std::cout << "\n\n" << std::endl;

  return 0;
}
