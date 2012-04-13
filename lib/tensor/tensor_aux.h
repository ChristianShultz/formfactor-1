#ifndef tensor_AUX_H_H_GUARD
#define tensor_AUX_H_H_GUARD

#include "tensor.h"

namespace tensor
{
  /**
     @file tensor_aux.h
     @brief defines useful functions for minimal #include 
   */

  /**
     @brief get a copy of the tensor from the base pointer

     @details Tensor derives from TensorImplbase which is pure virtual
     so we can use dynamic_cast to type check at runtime and break on 
     null/failure using assert(ptr).
   */
    template<typename T, idx_t rank>
    tensor::Tensor<T, rank> castAndDeref(const tensor::TensorImplBase *pVecPtr)
    {
        const tensor::Tensor<T, rank> *retPtr
            = dynamic_cast<const typename tensor::Tensor<T, rank>* >(pVecPtr);
        POW2_ASSERT(retPtr);                      // make sure the cast worked
        tensor::Tensor<T, rank> ret = *retPtr;    // copy the casted pointer into an actual object
        return ret;                               // return that object
    }

}
#endif
