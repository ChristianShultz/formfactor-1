#ifndef KINFAC_KEY_H_H_GUARD
#define KINFAC_KEY_H_H_GUARD

#include <string>
#include "utils/pow2assert.h"
#include "kinfac_fwd.h"

/**
   @file kinfac_key.h
   @brief the keys used for constructing the rows of the linear system
 */

namespace kinfac
{

    /**
       @brief a primitive key to hold src/sink information

     */
    struct KinKey_p
    {
        KinKey_p(void); //! not implemented
        KinKey_p(const EnsemReal &E_,
                 const Array<int> &mom_,
                 const idx_t J_,
                 const short helicity_,
                 const std::string &parity_)
            : E(E_) , mom(mom_) , J(J_) , helicity(helicity_) , parity(parity_)
        {
            POW2_ASSERT_DEBUG((parity == "+") || (parity == "-"));
            POW2_ASSERT_DEBUG(helicity <= J);
        }

        KinKey_p &operator=(const KinKey_p &o)
        {
            if(this != &o)
                {
                    E = o.E;
                    mom = o.mom;
                    J = o.J;
                    helicity = o.helicity;
                    parity = o.parity;
                }

            return *this;
        }

      KinKey_p(const KinKey_p &o)
      : E(o.E) , mom(o.mom) , J(o.J) , helicity(o.helicity) , parity(o.parity)
      {      }

        EnsemReal E;
        Array<int> mom;
        idx_t J;
        short helicity;
        std::string parity;
    };

    typedef KinKey_p SourceKey;
    typedef KinKey_p SinkKey;

    /**
       @brief a key to access kinematic factors
     */
    struct KinKey
    {
      KinKey(void); //! not implemented

      //! nb this is read right to left like a matrix element
      KinKey(const SinkKey &_sink, const SourceKey &_src) // now I'm coding right to left.. sigh
	: src(_src) , sink(_sink)
      { }

        SourceKey src;
        SinkKey sink;
    };


    /**
       @brief the registration isn't thread safe.
       @details we can pre-register the keys before parallelization
       so that we only have to do reads on the static map which should be
       thread safe, it is wrapped in an omp critical block
     */
    void registerPolarisationTensors(const KinKey &key);

} // close kinfac namespace

#endif
