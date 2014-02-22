/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : grab_ff.cc

* Purpose :

* Creation Date : 22-04-2013

* Last Modified : Sat 22 Feb 2014 05:05:47 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>

#include "radmat/ff/lorentzff_formfactor_factory.h"
#include "radmat/ff/lorentzff_canonical_rotations.h"
#include "radmat/utils/tensor.h"
#include "radmat/register_all/register_all.h"

//
//
// NB: this uses lorentzff_formfactor_factory
//
//  the internal form factors used by radmat are 
//  mapped to linear combinations of these but are
//  evaluated in a slightly different way, check 
//  radmat/ff_interface/forfactor_factory.h for
//  more details 
//
//
//  -- currently this guy does not have access to 
//  the cubic form factors 
//
//  -- currently this guy does not have access to 
//  the helicity form factors
//
//



using namespace radmat;
using namespace std; 

struct inp
{
  string matElemID; 
  Tensor<double,1> pf,pi;
  mom_t l,r; 
  double pfac;
  bool jmu;
  int mu, hf, hi; 
};


// assuming that the momentum factor is already in the space bits here -- on shell
void do_boost(Tensor<double,1> &p)
{
  p[0] = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
}


inp usage(int argc, char *argv[])
{
  if( (argc != 13) && (argc != 14) ) 
  {
    std::cerr << "usage: <grab_ff> <matElemID> <m_f> <nx ny nz>  <m_i> <nx ny, nz> <2pi/xi/nspace> <h_f> <h_i>  <optional jmu>" << std::endl;
    exit(1); 
  }


  inp foo;

  { std::stringstream ss(argv[1]); ss >> foo.matElemID; }
  { std::stringstream z(argv[2]),o(argv[3]),t(argv[4]),th(argv[5]);
    mom_t p; 
    p.resize(3); 
    Tensor<double,1> fred((TensorShape<1>())[4],0.);
    z >> fred[0];
    o >> fred[1];
    t >> fred[2];
    th >> fred[3];
    
    p[0] = fred[1]; 
    p[1] = fred[2];
    p[2] = fred[3];  

    foo.l = p; 
    foo.pf = fred; }

    { std::stringstream z(argv[6]),o(argv[7]),t(argv[8]),th(argv[9]);
    Tensor<double,1> fred((TensorShape<1>())[4],0.);
    mom_t p; 
    p.resize(3); 
    z >> fred[0];
    o >> fred[1];
    t >> fred[2];
    th >> fred[3];
     
    
    p[0] = fred[1]; 
    p[1] = fred[2];
    p[2] = fred[3];  

    foo.r = p; 
    foo.pi = fred; }

    {std::stringstream ss(argv[10]); ss >> foo.pfac;}
    {std::stringstream ss(argv[11]); ss >> foo.hf;}
    {std::stringstream ss(argv[12]); ss >> foo.hi;}


    for(int i = 1; i < 4; ++i)
    {
      foo.pf[i] *= foo.pfac;
      foo.pi[i] *= foo.pfac;
    }
    
    do_boost(foo.pf);
    do_boost(foo.pi); 

    foo.jmu = false; 

    if(argc == 14)
    {
      foo.jmu = true;
      std::stringstream ss(argv[11]); ss >> foo.mu; 
    }

    return foo; 
}


int main(int argc, char *argv[])
{

  AllFactoryEnv::registerAll(); 

  inp fred = usage(argc, argv); 

  std::cout << "the canonical frame is " <<
    radmat::LatticeRotationEnv::TheRotationGroupGenerator::Instance().get_can_frame_string(fred.l,fred.r) << std::endl;

  rHandle<FFAbsBase_t > foo = LorentzffFormFactorDecompositionFactoryEnv::callFactory(fred.matElemID); 

  itpp::Mat<std::complex<double> > baz = (*foo)( 
      std::pair< Tensor<double,1> , int >(fred.pf,fred.hf),
      std::pair< Tensor<double,1> , int >(fred.pi,fred.hi),
      fred.pfac); 

  if(fred.jmu)
  {

    itpp::Vec<std::complex<double> > bar;
    bar = itpp::round_to_zero(baz.get_row(fred.mu) , 1e-6); 

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

    std::cout << "by value (row index is lorentz, col is FF num) \n" 
      << itpp::round_to_zero( baz , 1e-6);


    std::cout << "\n\n" << std::endl;
  }

  return 0;
}
