#ifndef FAKE_H_H_GUARD
#define FAKE_H_H_GUARD

#include "../ff_base.h"
#include "radmat/utils/tensor.h"

namespace radmat
{
  namespace fake
  {
    struct ff_fake : public FFBlockBase_t
    {
      typedef FFBlockBase_t::pack_handle pack_handle;
   
      FFBlock_rt operator()(const pack_handle &final, const pack_handle &initial) const
      {
	Tensor<std::complex<double> , 1 > pp, pm;
	pp = final->pmu() + initial->pmu();   // test the asymmetric assignment
	pm = final->pmu() - initial->pmu();
	std::complex<double> Q2 = -pm*(g_dd()*pm);
	std::complex<double> mm2 = pp*(g_dd()*pm);


	return ( (Q2/mm2)*pp + pm);
      }
    };
  }
}
#endif
