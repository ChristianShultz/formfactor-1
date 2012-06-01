#ifndef FF_BASE_H_H_GUARD
#define FF_BASE_H_H_GUARD

#include "xml_array.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/pow2assert.h"
#include "ensem/ensem.h"
#include "adat/handle.h"
#include <utility>


/**
   @file ff_base.h
   @brief some common definitions and pure virtual abstract base classes
 */


namespace radmat
{

  /**
     @brief hold the lorentz invariants (info to make them) of one side of a matrix element
     @details we will used shared pointers to these to allow easy updating of the only cfg dependent
              paramater, the energy of the meson.  This allows us to repeatedly use the single cfg 
	      linear system row generator to spit out the info for all cfgs via operator()(double Energy_cfg).
   */

    struct InvariantPack_t
    {
        InvariantPack_t(const double _E, const XMLArray::Array<int> &_mom, const idx_t _J, const short _lambda)
            : E(_E) , mom(_mom) , J(_J), lambda(_lambda)
        {
            POW2_ASSERT(abs(lambda) <= J);
        }

        void operator()(const double _E) // update the energy for looping
        {
            E = _E;
        }

        Tensor<double, 1> pmu(void)
        {
            Tensor<double, 1> p(TensorShape<1>()[4], 0.);
            p[0] = E;
            p[1] = mom[0];
            p[2] = mom[1];
            p[3] = mom[2];
            return p;
        }

        double E;
        XMLArray::Array<int> mom;
        idx_t J;
        short lambda;
    };

  // since everything except <0|j^{\mu}|0> involves polarisation tensors we will only 
  // bother to construct complex linear systems 
    typedef Tensor<std::complex<double>, 1 > FFBlock_rt;

  //! a pure virtual abc for a single term in the form factor decomposition, a list of these forms a decomp
    struct FFBlockBase_t
    {
        typedef ADAT::Handle<InvariantPack_t> pack_handle;
        virtual FFBlock_rt operator()(const pack_handle &final, const pack_handle &initial) const = 0;
    };

}// namespace


#endif
