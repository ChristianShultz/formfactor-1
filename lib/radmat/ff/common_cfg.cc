// common_cfg.cc -
//
// Wednesday, May 30 2012
//

#include"common_cfg.h"

namespace radmat
{

  Tensor<double,1> pPlus(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i)
  {
    return p_f + p_i;
  }

  Tensor<double,1> pMinus(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i)
  {
    return p_f - p_i;
  }

} // close radmat
