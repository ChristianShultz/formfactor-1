// test_clebsch.cc -
//
// Thursday, March 29 2012
//

#include "hadron/clebsch.h"
#include "utils/pow2assert.h"
#include <iostream>

using namespace Hadron;

int 
main(void)
{
  POW2_ASSERT(clebsch(1,1,1,1,2,2) == 1.);
  std::cout << "Adat clebsch successfully linked and tested" << std::endl;
}
