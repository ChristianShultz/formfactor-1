// tensor_algebra.cc -
//
// Tuesday, April 10 2012
//

#include <vector>
#include"tensor_algebra.h"

using namespace tensor;

Tensor<short,2> tensor::g_uv(void)
{
  Tensor<short,2> ret;
  ret.create(std::vector<tensor::idx_t>(2,4));

  ret[0][0] = 1;
  ret[1][1] = -1;
  ret[2][2] = -1;
  ret[3][3] = -1;
  
  return ret;
}
