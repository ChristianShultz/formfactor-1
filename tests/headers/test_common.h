#ifndef TEST_COMMON_H_H_GUARD
#define TEST_COMMON_H_H_GUARD

#include "tester.h"

namespace radmat
{

  tester test_polarisation_tensor(void);
  tester test_ff_fake(void);

  namespace PiPi
  {
    tester test_ff_debug(void);
    tester test_ff(void);
  }

}
#endif
