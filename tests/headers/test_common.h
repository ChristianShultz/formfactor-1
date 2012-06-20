#ifndef TEST_COMMON_H_H_GUARD
#define TEST_COMMON_H_H_GUARD

#include "tester.h"
#include <string>
#include <complex>

namespace radmat
{
  // utils
  tester test_polarisation_tensor(void);
  tester test_ff_fake(void);

  // formfacs
  namespace PiPi
  {
    tester test_ff_debug(void);
    tester test_ff(void);
  }

  // llsq 
  tester test_LLSQ_solver_SVDMakeSquare(void);

  tester test_solver_factory(void);
  tester test_mat_elem_factory(void);

  // fake data
  struct test_covarrying_vectors
  {
    template<typename T>
    tester test(void) const;
  };

  template<> tester test_covarrying_vectors::test<double>(void) const;
  template<> tester test_covarrying_vectors::test<std::complex<double> >(void) const;

  tester test_minimal_fake_data(const std::string &matElemID);

  

}
#endif
