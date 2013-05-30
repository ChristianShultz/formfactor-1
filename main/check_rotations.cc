/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : check_rotations.cc

 * Purpose :

 * Creation Date : 07-05-2013

 * Last Modified : Fri May 10 09:05:31 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/utils/tensor.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/polarisation_tensors.h"
#include "hadron/irrep_util.h"
#include <vector>
#include <list>
#include <iostream>
#include <complex>

using namespace radmat; 


  Tensor<std::complex<double> , 1> 
toTens(const XMLArray::Array<int> &p, const double E, const double pfac)
{
  std::complex<double> one(1.,0.); 
  Tensor<std::complex<double> , 1> ret((TensorShape<1>())[4],0.);
  ret[0] = E*one; 
  ret[1] = (pfac*p[0])*one; 
  ret[2] = (pfac*p[1])*one; 
  ret[3] = (pfac*p[2])*one; 

  return ret; 
}

  Tensor<double, 1> 
toTens(const XMLArray::Array<int> &p)
{
  Tensor<double,1> ret((TensorShape<1>())[3],0.); 
  ret[0] = double(p[0]);
  ret[1] = double(p[1]);
  ret[2] = double(p[2]);

  return ret; 
}

  Tensor<double, 1> 
zTens(const XMLArray::Array<int> &p)
{
  Tensor<double,1> ret((TensorShape<1>())[3],0.); 
  double pp =  double(p[0])* double(p[0]) + double(p[1])* double(p[1]) + double(p[2])* double(p[2]);
  ret[2] = sqrt(pp); 
  return ret; 
}


  void 
print (const XMLArray::Array<int> &mom)
{
  Tensor<double,1> p = toTens(mom); 
  Tensor<double,1> pz = zTens(mom); 

# if 0

  std::cout << p << "   " 
    << pz << "   "
    << genRotationMatrix3D(mom) * pz << "   "
    << p - genRotationMatrix3D(mom) * pz  << std::endl;

#endif 

  Tensor<double,1> flag = p - genRotationMatrix3D(mom) * pz ; 

  for(int i = 0; i < 3; ++i)
  {
    POW2_ASSERT(fabs(flag[i]) < 0.00001);  
  }



  double E = 0.2161; // 743 rho mass
  double pfac = 0.11; // 2pi/xi/Ls
  int J = 1; 

  Tensor<std::complex<double> , 1> pmu = toTens(mom,E,pfac); 
  genPolTens<1> eps_gen(mom); 

  // sanity check on polarization tensors -- must be orthogonal 
  for(short lambda = -J; lambda <= J; ++lambda)
  {
    std::complex<double> foo = (pmu * ( g_dd() * eps_gen(E,lambda,pfac) ) ).value() ; 

    POW2_ASSERT(std::norm(foo) < 0.00001);  
  }

}


  itpp::Mat<std::complex<double> > 
eps3d(const ADATXML::Array<int> &mom , const bool create)
{
  Tensor<std::complex<double>, 1 > tmp;
  genPolTens3D<1> eps(mom);
  itpp::Mat<std::complex<double> > eps3(3,3); 

  for(int h = 1; h > -2; --h)
  {
    tmp = eps.get(h);

    for(int i = 0; i < 3; ++i)
      if(create)
        eps3(1-h,i) = tmp[i];
      else
        eps3(1-h,i) = std::conj(tmp[i]); 

  }

  return eps3; 
}

void
print3d(const XMLArray::Array<int> &mom)
{
  std::cout << "p = " << mom[0] << mom[1] << mom[2] << std::endl;

  std::cout << itpp::round_to_zero(eps3d(mom,true),0.00001) << std::endl;
}



  int 
main(void)

{

  std::list<XMLArray::Array<int> > D4,D2,D3;  

  { XMLArray::Array<int> foo; foo.resize(3); foo[0] = 1; foo[1] = 0; foo[2] = 0; D4 = Hadron::generateLittleGroupMom("D4",foo); }
  { XMLArray::Array<int> foo; foo.resize(3); foo[0] = 1; foo[1] = 1; foo[2] = 0; D2 = Hadron::generateLittleGroupMom("D2",foo); }
  { XMLArray::Array<int> foo; foo.resize(3); foo[0] = 1; foo[1] = 1; foo[2] = 1; D3 = Hadron::generateLittleGroupMom("D3",foo); }


  std::list<XMLArray::Array<int> >::const_iterator it; 

  std::cout << " checking the following .. will crash if something fails" << std::endl;

  std::cout << "********* D4 " << std::endl;
  for(it = D4.begin(); it != D4.end(); ++it)
    print(*it); 

  std::cout << "********* D2 " << std::endl;
  for(it = D2.begin(); it != D2.end(); ++it)
    print(*it); 



  std::cout << "********* D3 " << std::endl;
  for(it = D3.begin(); it != D3.end(); ++it)
    print(*it); 



  std::cout << "********* Checking 3d " << std::endl;


  for(it = D4.begin(); it != D4.end(); ++it)
    print3d(*it); 


  return 0; 
}
