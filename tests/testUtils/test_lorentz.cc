// test_lorentz.cc -
//
// Thursday, April 19 2012
//

#include "utils/lorentz.h"
#include "utils/pow2assert.h"
#include "itpp/itbase.h"

int
main(void)
{

  double precision = 1e-10; 

  double m2 = 3.14159265;         
  itpp::Vec<double> pmu(4);
  pmu.zeros();

  pmu[0] = sqrt(m2);  

  double pz_boost = 1;
  double E_boost = sqrt(m2 + pz_boost*pz_boost); 
  double rapidity_boost = lorentz::rapidity(E_boost, pz_boost); 

  itpp::Mat<double> boost_z = lorentz::boostZ(rapidity_boost);
  itpp::Vec<double> pmu_boostz = boost_z*pmu;
  
  /*
  std::cout << boost_z << std::endl;
  std::cout << pmu << std::endl;
  std::cout << pmu_boostz << std::endl;
  */
  
  // check that we are still on the mass shell
  POW2_ASSERT( fabs(pmu_boostz*(lorentz::gmunu()*pmu_boostz) - m2) < precision);
 
  itpp::Mat<double> rotate_to_minus_z = lorentz::rotateToZAxis(pmu_boostz, false);
  itpp::Mat<double> boost_minus_z = rotate_to_minus_z*boost_z;
  itpp::Vec<double> pmu_boost_minus_z = boost_minus_z*pmu; 

  // check that the rotation left us on mass shell
  POW2_ASSERT( fabs(pmu_boost_minus_z*(lorentz::gmunu()*pmu_boost_minus_z) - m2)  < precision ); 

  /*
  std::cout << boost_z*pmu << std::endl; 
  std::cout << boost_minus_z*pmu << std::endl; 
  */
  
  itpp::Vec<double> p(4), p_prime(4);
  
  double m = sqrt(m2);
  double m_prime = m2;                       
  
  p = itpp::randu(4); 
  p_prime = itpp::randu(4);  

  std::cout << "m_p " << m << std::endl;
  std::cout << "m_p_prime " << m_prime << std::endl;

  std::cout << "p " << p << std::endl; 
  std::cout << "p_prime " << p_prime << std::endl; 

  p[0] = sqrt(m*m + p[1]*p[1] + p[2]*p[2] + p[3]*p[3] ); 
  p_prime[0] = sqrt(m_prime*m_prime + p_prime[1]*p_prime[1] 
		    +p_prime[2]*p_prime[2] + p_prime[3]*p_prime[3] ); 


  lorentz::LorentzTransform LT = lorentz::genBreitLT(p_prime, p);   

   
  return 0;
}
