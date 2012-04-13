#ifndef tensor_ALGEBRA_H_H_GUARD
#define tensor_ALGEBRA_H_H_GUARD


#include "tensor.h"
#include "narray.h"
#include <iostream>
#include <vector>
#include "utils/pow2assert.h"
#include "ensem/ensem.h"


namespace tensor
{
  /**
     @file tensor_algebra.h
     @brief contains templated algebraic operations for tensors
   */

  //! \f$ g^{\mu\nu} \f$
  Tensor<short,2> g_uv(void);

  /** \addtogroup tensor_algebra
   * @{
   */

    // fwd

    //! overload * operator for rank 1 tensors to make it return a scalar
    template<typename T, typename U>
    typename ENSEM::Promote<T, U>::Type_t operator*(const Tensor<T, 1> &lhs, const Tensor<U, 1> &rhs);

    //! overload contract for rank 1 tensors to return a scalar since rank 0 is ill defined
    template<typename T, typename U>
    typename ENSEM::Promote<T, U>::Type_t contract(const Tensor<T, 1> &lhs, const Tensor<U, 1> &rhs);

  /** 
      @brief overload * operator for arbitrary rank tensors to contract the last index of lhs with the
      first index of rhs
  */
    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M - 2 > operator*(const Tensor<T, N> &lhs, const Tensor<U, M> &rhs);

    //! tensor product to build up bigger rank tensors from smaller ones
    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M > tensorProduct(const Tensor<T, N> &A, const Tensor<U, M> &B);

    //! overload wedge for tensor product
    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M > operator^(const Tensor<T, N> &a, const Tensor<U, M> &b);

    //! overload ==
    template<typename T, typename U, idx_t N>
    bool operator==(const Tensor<T, N> &a, const Tensor<U, N> &b);

    //! overload !=
    template<typename T, typename U, idx_t N>
    bool operator!=(const Tensor<T, N> &a, const Tensor<U, N> &b);


  /** @} */ // tensoralgebra group


  /** \addtogroup levi_civita
   * @{
   */

    //! determine if elements of ptr are unique
    template<idx_t N>
    bool filter(const  idx_t *const ptr);

    //! bubble sort an index a,b,c..n to determine the number of permutations from  0,1,2..n 
    template<idx_t N>
    idx_t leviBubbleSort(const idx_t *const ptr);

    //! generate an n-rank levi civita tensor
    template<typename T, idx_t N>
    Tensor<T, N> leviCivitaTensor(void);

  /** @} */ // levicivita group

    // impl
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////

    template<typename T, typename U>
    typename ENSEM::Promote<T, U>::Type_t operator*(const Tensor<T, 1> &lhs, const Tensor<U, 1> &rhs)
    {
        idx_t dim = lhs.getDim(0);
        pow2assert(dim == rhs.getDim(0));
        pow2assert(lhs.isRaised(0) != rhs.isRaised(0));

        typename ENSEM::Promote<T, U>::Type_t ret;
        ret = 0.;

        for(idx_t i = 0; i < dim; ++i)
            ret += lhs[i] * rhs[i];

        return ret;
    }

    template<typename T, typename U>
    typename ENSEM::Promote<T, U>::Type_t contract(const Tensor<T, 1> &lhs, const Tensor<U, 1> &rhs)
    {
        return lhs * rhs;
    }

    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M - 2 > operator*(const Tensor<T, N> &lhs, const Tensor<U, M> &rhs)
    {
        return contract(lhs, rhs, N - 1, 0);
    }

    // this is a really really dumb way to store the tensor product but it fits easily into the Tensor framework
    // and we likely don't need to go beyond polarization tensors for spin-4 so it isnt a big deal
    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M > tensorProduct(const Tensor<T, N> &A, const Tensor<U, M> &B)
    {
        Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M > ret;
        std::vector<idx_t> dim;

        for(idx_t i = 0; i < N; ++i)
            dim.push_back(A.getDim(i));

        for(idx_t i = 0; i < M; ++i)
            dim.push_back(B.getDim(i));

        POW2_ASSERT(ret.create(&dim[0]));

        idx_t **a, **b, **r, *cycle;

        a = new idx_t*[N];
        b = new idx_t*[M];
        r = new idx_t*[N + M];
        cycle = new idx_t[N + M];

        for(idx_t i = 0; i < N; ++i)
            {
                cycle[i] = 0;
                a[i] = &cycle[i];
                r[i] = &cycle[i];
            }

        for(idx_t i = 0; i < M; ++i)
            {
                cycle[N + i] = 0;
                b[i] = &cycle[N + i];
                r[N + i] = &cycle[N + i];
            }

        // standard n-dimensional for loop construction
        while(true)
            {

                ret.getElem(r) = A.getElem(a) * B.getElem(b);

                idx_t j;

                for(j = 0; j < N + M; ++j)
                    {
                        ++cycle[j];

                        if(cycle[j] < dim[j])
                            break;

                        cycle[j] = 0;
                    }

                if(j == N + M) // break while
                    break;
            }

        delete [] a;
        delete [] b;
        delete [] cycle;

        return ret;
    }

    // overload wedge for tensor product
    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t, N + M > operator^(const Tensor<T, N> &a, const Tensor<U, M> &b)
    {
        return tensorProduct(a, b);
    }

    // ==
    template<typename T, typename U, idx_t N>
    bool operator==(const Tensor<T, N> &a, const Tensor<U, N> &b)
    {
        idx_t *dim = new idx_t[N];
        idx_t *idx = new idx_t[N];
        idx_t **pidx = new idx_t*[N];

        for(idx_t i = 0; i < N; ++i)
            {
                dim[i] = a.getDim(i);

                if(dim[i] != b.getDim(i))
                    return false;

                idx[i] = 0;
                pidx[i] = &idx[i];
            }

        while(true)
            {
                if(a.getElem(pidx) != b.getElem(pidx))
                    return false;

                idx_t j;

                for(j = 0; j < N; ++j)
                    {
                        ++idx[j];

                        if(idx[j] < dim[j])
                            break;

                        idx[j] = 0;
                    }

                if(j == N)
                    break;
            }

        delete [] dim;
        delete [] idx;
        delete [] pidx;

        return true;
    }


    template<typename T, typename U, idx_t N>
    bool operator!=(const Tensor<T, N> &a, const Tensor<U, N> &b)
    {
        return !(a == b);
    }

    // levi civita helpers

    // filter out the thing like eps_iijk
    template<idx_t N>
    bool filter(const  idx_t *const ptr)
    {
        std::vector<idx_t> ct(N, 0);

        for(idx_t i = 0; i < N; ++i)
            {
                if(ct[ptr[i]] == 0)
                    ++ct[ptr[i]];
                else
                    return false;
            }

        return true;
    }

    // this will count the number of permutations from 0,1,2,3,..,n
    // there must be a smarter way but this is incredibly simple and
    // easy to implement, the tensor element is then just -1^nperm
    template<idx_t N>
    idx_t leviBubbleSort(const idx_t *const ptr)
    {
        idx_t *idxs = new idx_t[N];
        idx_t nperm = 0, tmp;

        memcpy(idxs, ptr, N * sizeof(idx_t));

        for(idx_t i = 0; i < N; ++i)                 // count the number of permutations to put
            for(idx_t j = 0; j < i; ++j)               // the set of indicies in the order
                if(idxs[i] < idxs[j])                    // 0,1,2,..n
                    {
                        ++nperm;
                        tmp = idxs[i];
                        idxs[i] = idxs[j];
                        idxs[j] = tmp;
                    }

        delete [] idxs;

        return nperm;
    }

    template<typename T, idx_t N>
    Tensor<T, N> leviCivitaTensor(void)
    {
        Tensor<T, N> levi;
        std::vector<idx_t> dim(N, N);
        POW2_ASSERT(levi.create(&dim[0]));

        idx_t *idx = new idx_t[N];

        for(idx_t i = 0; i < N; ++i)
            idx[i] = 0;

        // the standard n-dimensional for loop we've been using
        while(true)
            {
                if(!!!filter<N>(idx))  // find the zero elements 
                    levi.getElem(idx) = T(0.);
                else
                    {
                        if((leviBubbleSort<N>(idx) % 2) == 0)
                            levi.getElem(idx) = T(1);
                        else
                            levi.getElem(idx) = T(-1);

                        /*  DEBUG

                          std::cout << "nperm = " << leviBubbleSort<N>(idx) <<  " | nperm % 2 = " << leviBubbleSort<N>(idx) % 2 << " | ";

                          for(int i = 0; i < N; ++i)
                          std::cout << idx[i] << " ";

                          std::cout << " " <<  levi.getElem(idx) << std::endl;

                        */

                    }

                // cycle through indicies -- leftmost index cycles fastest
                idx_t j;

                for(j = 0; j < N; ++j)
                    {
                        ++idx[j];

                        if(idx[j] < N)
                            break;

                        idx[j] = 0;
                    }

                // break while loop
                if(j == N)
                    break;
            }

        delete [] idx;

        return levi;
    }

} // close tensor namespace

#endif
