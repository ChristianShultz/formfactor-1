#ifndef tensor_FWD_H_H_GUARD
#define tensor_FWD_H_H_GUARD

namespace tensor
{
  /**
     @file tensor_fwd.h
     @brief forward declaration of some structs
   */
    typedef unsigned short idx_t;

    template<typename T, idx_t N>
    struct sub_t;

    template<typename T, idx_t N>
    struct narray_t;

    // unary minus
    template<typename T, idx_t N>
    narray_t<T, N> operator-(const narray_t<T, N> &plus);

    template<typename T, idx_t N>
    struct Tensor;

    // unary minus
    template<typename T, idx_t N>
    Tensor<T, N> operator-(const Tensor<T, N> &plus);
}
#endif
