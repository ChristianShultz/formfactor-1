// common_cfg.cc -
//
// Wednesday, May 30 2012
//

#include"common_cfg.h"

namespace radmat
{

  Tensor<double,1> pPlus(const PInv_t &moms)
  {
    return moms.first + moms.second;
  }

  Tensor<double,1> pMinus(const PInv_t &moms)
  {
    return moms.first - moms.second;
  }

} // close radmat
