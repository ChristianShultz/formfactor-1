/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : grab_ff.cc

* Purpose :

* Creation Date : 22-04-2013

* Last Modified : Tue 15 Oct 2013 07:48:41 PM EDT

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
  bool jmu;
  int mu; 
};


// assuming that the momentum factor is already in the space bits here -- on shell
void do_boost(Tensor<double,1> &p)
{
  p[0] = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
}


inp usage(int argc, char *argv[])
{
  if( (argc != 11) && (argc != 12) ) 
  {
    std::cerr << "usage: <grab_ff> <matElemID> <m_f> <nx ny nz>  <m_i> <nx ny, nz> <2pi/xi/nspace>  <optional jmu>" << std::endl;
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
    
    do_boost(foo.pf);
    do_boost(foo.pi); 

    foo.jmu = false; 

    if(argc == 12)
    {
      foo.jmu = true;
      std::stringstream ss(argv[11]); ss >> foo.mu; 
    }

    return foo; 
}


int main(int argc, char *argv[])
{

  inp fred = usage(argc, argv); 

  FormFactorDecompositionFactoryEnv::registerAll(); 

  rHandle<ffBase_t<std::complex<double> > > foo = FormFactorDecompositionFactoryEnv::callFactory(fred.matElemID); 

  if(fred.jmu)
  {

    itpp::Mat<std::complex<double> > baz = (*foo)(fred.pf,fred.pi,fred.pfac); 
    itpp::Vec<std::complex<double> > bar = baz.get_row(fred.mu); 

    for(int elem = 0; elem < bar.size() -1; ++elem)
    {
      std::complex<double> my_cmp = bar[elem];
      std::cout << my_cmp.real() << " " << my_cmp.imag() << "   ,   "; 
    }

    std::complex<double> my_cmp = bar[bar.size() -1];
    std::cout << my_cmp.real() << " " << my_cmp.imag() << std::endl;

  }
  else
  {
    std::cout << " The list of ffs is : \n" << foo->ff() << std::endl;

    std::cout << "by value (row index is lorentz, col is FF num) \n" << (*foo)(fred.pf,fred.pi,fred.pfac);

    std::cout << "\n\n" << std::endl;
  }

  return 0;
}
