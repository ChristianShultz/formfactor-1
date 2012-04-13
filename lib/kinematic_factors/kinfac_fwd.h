#ifndef KINFAC_FWD_H_H_GUARD
#define KINFAC_FWD_H_H_GUARD


#include "ensem/ensem.h"
#include "xml_array.h"
#include "tensor/tensor_fwd.h"
#include "itpp/itbase.h"
#include "tensor/tensor.h"

#include <vector>

/**
   @file kinfac_fwd.h
   @brief hold some useful fwd information
 */

namespace kinfac
{

    // bring EnsemReal/Complex into this namespace
    using ENSEM::EnsemReal;
    using ENSEM::EnsemComplex;
    using XMLArray::Array;
    using tensor::idx_t;

}

#endif
