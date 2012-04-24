// lorentz.cc -
//
// Thursday, April 19 2012
//


#include "lorentz.h"
#include <math.h>
#include "pow2assert.h"
 

using namespace lorentz;

double lorentz::rapidity(const double E, const double pz)
{
  return 0.5*log( (E - pz)/(E + pz) );
}

itpp::Mat<double> lorentz::rotateToZAxis(const itpp::Vec<double> &p4, const bool plus)
{
  if(plus) 
    return rotMatrix4( genEulerAnglesToZAxis(p4) );

  return rotMatrix4( genEulerAnglesToNegativeZAxis(p4) ); 
}

itpp::Mat<double> lorentz::boostZ(const double xi /*rapidity*/)
{
  itpp::Mat<double> K_3(4,4);
  K_3.zeros();
  itpp::Vec<double> I_4(4);
  I_4.ones();

  K_3(3,0) = 1;
  K_3(0,3) = 1;

  return  (
	   itpp::diag(I_4) 
	   - K_3*K_3 
	   - K_3*sinh(xi) 
	   + K_3*K_3*cosh(xi)
	   );    
}


/*
  x-convention. 
  In this convention the rotation is given by Euler angles (a,b,c) , where

  1. the first rotation is by an angle 'c' about the z-axis using C,

  2. the second rotation is by an angle 'b' about the former x-axis using B, and

  3. the third rotation is by an angle 'a' about the former z-axis using A.

  the return is a lorentz rotation matrix with the O(3) matrix embedded in the spatial section 
  and a 1 in the 0,0 temporal element
*/

itpp::Mat<double> lorentz::rotMatrix4(const itpp::Vec<double> &euler_angles)
{
  itpp::Mat<double> A(4,4) , B, C;
  A.zeros(); 
  A(0,0) = 1;
  A(1,1) = 1;
  A(2,2) = 1;
  A(3,3) = 1;
  B = A; 
  C = A;

  A(1,1) = cos(euler_angles[0]); 
  A(2,2) = cos(euler_angles[0]);
  A(1,2) = sin(euler_angles[0]);
  A(2,1) = -sin(euler_angles[0]);
  
  B(2,2) = cos(euler_angles[1]); 
  B(3,3) = cos(euler_angles[1]); 
  B(2,3) = sin(euler_angles[1]); 
  B(3,2) = -sin(euler_angles[1]); 

  C(1,1) = cos(euler_angles[2]); 
  C(2,2) = cos(euler_angles[2]);
  C(1,2) = sin(euler_angles[2]);
  C(2,1) = -sin(euler_angles[2]);

  // std::cout << "A" << A << std::endl;
  // std::cout << "B" << B << std::endl;
  // std::cout << "C" << C << std::endl;

  return A*B*C;
}
  
itpp::Vec<double> lorentz::genEulerAnglesToZAxis(const itpp::Vec<double> &p4)
{

  POW2_ASSERT(p4.size() == 4);

  itpp::Vec<double>  eul(3);

  // need to be careful to keep the inputs to atrig fcns representable
  
  const static double pi = acos(-1.0);

  // last rotation wouldn't do anything since its on the z axis
  // be careful of phases later for the wigned D matrix elems
  eul[0] = 0;

  double r_y = sqrt(p4[1]*p4[1] + p4[2]*p4[2]);

  if(p4[3] == 0.)
    eul[1] = 0.5*pi;
  else
    eul[1] = p4[3] < 0. ? pi -atan(r_y/p4[3]) : -atan(r_y/p4[3]); //r lies in the positive y direction
  
  if(p4[2] == 0.)
    eul[2] = p4[1] < 0. ? -0.5*pi : 0.5*pi; 
  else
    eul[2] = p4[2] < 0. ? pi -atan(p4[1]/p4[2]) : -atan(p4[1]/p4[2]);
      
  return eul;
}

itpp::Vec<double> lorentz::genEulerAnglesToNegativeZAxis(const itpp::Vec<double> &v4)
{
  itpp::Vec<double> eul = genEulerAnglesToZAxis(v4); 
  eul[1] += acos(-1.0); // rotate by compliment angle
  return eul; 
}


itpp::Mat<double> lorentz::gmunu(void) 
{
  itpp::Vec<double> g(4);
  g.ones(); 
  g = -g;
  g[0] = 1;
  return itpp::diag(g); 
}

lorentz::LorentzTransform lorentz::genBreitLT(const itpp::Vec<double> &pp4, const itpp::Vec<double> &p4)
{
  itpp::Mat<double> Rp4z = lorentz::rotateToZAxis(p4); 
  itpp::Vec<double> p4z = Rp4z*p4;
  double p4ZBoostRapidity = -lorentz::rapidity(p4z[0],p4z[3]);
  itpp::Mat<double> BzRp4z = lorentz::boostZ(p4ZBoostRapidity); 
  itpp::Mat<double> RppBzRp = lorentz::rotateToZAxis(BzRp4z*Rp4z*pp4,false); 
  itpp::Vec<double> RppBzRp_pp4 = RppBzRp*BzRp4z*Rp4z*pp4;
  double m = sqrt(p4[0]*p4[0] - p4[1]*p4[1] - p4[2]*p4[2] - p4[3]*p4[3]); 
  double E = RppBzRp_pp4[0];
  double pz = RppBzRp_pp4[3]; 
  double x = (E+m)/pz;
  double breitRapidity = 0.5* log( (x + 1.) / (x - 1.) ); // inverse hyperbolic cotangent

  lorentz::LorentzTransform foobar;
  foobar.Lambda = lorentz::boostZ(breitRapidity) * RppBzRp * BzRp4z * Rp4z;
  foobar.LambdaInv = (Rp4z.T()) * lorentz::boostZ(-p4ZBoostRapidity) * (RppBzRp.T()) * lorentz::boostZ(-breitRapidity); 

  return foobar; 
}
