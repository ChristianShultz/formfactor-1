// test_Levi.cc -
//
// Friday, March 30 2012
//

#include "radmat/tensor/tensorbase.h"
#include "radmat/utils/pow2assert.h"
#include <iostream>

using namespace radmat;

int factorial(int n)
{
  return (n == 1 ? 1 : n*factorial(n-1));
}


int 
main(void)
{

  int s2(0),s3(0),s4(0);

  Tensor<int,2> levi2 = leviCivitaTensor<int,2>();

  for(idx_t i = 0; i < 2; ++i)
    for(idx_t j = 0; j < 2; ++j)
      {
	POW2_ASSERT(levi2[i][j] == (j-i));
	s2 += levi2[i][j]*levi2[i][j];
      }
  POW2_ASSERT(s2 == factorial(2));


  Tensor<int,3> levi3 = leviCivitaTensor<int,3>();

  for(idx_t i = 0; i < 3; ++i)
    for(idx_t j = 0; j < 3; ++j)
      for(idx_t k = 0; k < 3; ++k)
	{
	  POW2_ASSERT(levi3[i][j][k] == (j-i)*(k-i)*(k-j)/2);
	  s3 += levi3[i][j][k]*levi3[i][j][k];
	}
  POW2_ASSERT(s3 == factorial(3));

  Tensor<int,4> levi4 = leviCivitaTensor<int,4>();
  for(idx_t i = 0; i < 4; ++i)
    for(idx_t j = 0; j < 4; ++j)
      for(idx_t k = 0; k < 4; ++k)
	for(idx_t l = 0; l < 4; ++l)
	  {
	    POW2_ASSERT(levi4[i][j][k][l] == (j-i)*(k-i)*(l-i)*(k-j)*(l-j)*(l-k)/12);
	    s4 += levi4[i][j][k][l]*levi4[i][j][k][l];
	  }
  POW2_ASSERT(s4 == factorial(4));

  std::cout << "Levi-Civita Tensor tested successfully" << std::endl;

  return 0;
}
