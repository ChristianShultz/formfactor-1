// test_INV.cc -
//
// Monday, April  2 2012
//

#include "radmat/polarisation/polarisationbase.h"
#include "radmat/tensor/tensorbase.h"
#include "radmat/utils/breit_frame.h"

using namespace radmat;
using namespace radmat::breit;

int
main(void)
{
  short h = -1;
  double precision = 1e-13;
  Tensor<double,1> pmu, p_z;
  std::vector<idx_t> dim(1,4);
  pmu.create(&dim[0]);

  pmu[1] = itpp::randu() + itpp::randu();
  pmu[2] = itpp::randu() + itpp::randu();
  pmu[3] = itpp::randu() + itpp::randu();


  p_z = pmu;
  p_z[1] = 0.;
  p_z[2] = 0.;
  p_z[3] = sqrt(pmu[1]*pmu[1] + 
		pmu[2]*pmu[2] + 
		pmu[3]*pmu[3]
		);

  itpp::Mat<double> R = rodRotMat(p_z,pmu);

  for(idx_t i = 0; i < 3; ++i)
    {
      POW2_ASSERT( fabs(pmu[i+1] - R(i,2)*p_z[3]) < precision);
    }
  pFacKey k(pmu,1,h), k_z(p_z,1,h);
  pFac factory;

  Coupler<1> c(k,&factory), c_z(k_z,&factory);
      
  Tensor<std::complex<double> , 1> *eps , *eps_z;
      
  // need to be called in this order to avoid a leak with valgrind b/c of 
  // operator() definition in Coupler<1>
  eps_z = c_z();
  eps = c();

  pFacInv::pFacInvRegister(k,eps);
  
  pFacInv::dumpInventory();

  std::cout << "results of test_INV in pFacInvDump.log " << std::endl;

  return 0;
}
