/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : grab_pol.cc

 * Purpose :

 * Creation Date : 22-04-2013

 * Last Modified : Mon Apr 22 16:29:44 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/ff/lorentzff_polarization_embedding.h"
#include <iostream>
#include <sstream>

using namespace radmat;
using namespace std;

void printTensor(const Tensor<double,1> &p, const int J, const int hel, const double pfac)
{
  switch(J)
  {
    case 1:
      {
      HelicityPolarizationTensor<1> foo;
      std::cout << foo(p,hel,pfac);
      break;
      }
    case 2:
      {
      HelicityPolarizationTensor<1> foo;
      std::cout << foo(p,hel,pfac);
      break;
      }
    case 3:
      {
      HelicityPolarizationTensor<2> foo;
      std::cout << foo(p,hel,pfac);
      break;
      }
    case 4:
      {
      HelicityPolarizationTensor<3> foo;
      std::cout << foo(p,hel,pfac);
      break;
      }
    default: 
      std::cout << "J = " << J << " is not supported" << std::endl;
      exit(1); 
  }
}


struct inp
{
  int J;
  int H; 
  Tensor<double,1> p;
  double pfac;
};



inp usage(int argc, char *argv[])
{
  if(argc != 8)
  {
    std::cerr << "usage: grab_pol <J> <H> <E> <nx ny nz> <2pi/xi/ns>" << std::endl;
    exit(1);
  }

  inp foo;

  {std::stringstream ss(argv[1]); ss >> foo.J;}

  {std::stringstream ss(argv[2]); ss >> foo.H;}

  {std::stringstream ss(argv[7]); ss >> foo.pfac;}

    { std::stringstream z(argv[3]),o(argv[4]),t(argv[5]),th(argv[6]);
    Tensor<double,1> fred((TensorShape<1>())[4],0.);
    z >> fred[0];
    o >> fred[1];
    t >> fred[2];
    th >> fred[3];
     
    foo.p = fred; }

    for(int i = 1; i < 4; ++i)
      foo.p[i] *= foo.pfac;

    return foo;
}



int main(int argc, char *argv[])
{

  inp foo = usage(argc, argv); 

  printTensor(foo.p,foo.J,foo.H,foo.pfac); 

  return 0; 
}
