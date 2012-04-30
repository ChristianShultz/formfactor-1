#ifndef NARRAY_STREAM_H_H_GUARD
#define NARRAY_STREAM_H_H_GUARD

#include "narray.h"
#include <iostream>

namespace radmat
{
  /**
     @file narray_stream.h

     @brief specializes printing for narray_t class
   */

  /**
     @brief A container for a function to allow us to specialize printing 
   */
    template<typename T, idx_t N>
    struct narray_stream
    {
        narray_stream(void) {}

      //! default implementation for N > 2, gets specialized for N = 1 or 2
        std::ostream &operator()(std::ostream &o, const narray_t<T, N> &arr);
    };


    template<typename T, idx_t N>
    std::ostream &narray_stream<T, N>::operator()(std::ostream &o, const narray_t<T, N> &a)
    {
        idx_t *ptr = new idx_t[N];

        for(idx_t i = 0; i < N; ++i)
            ptr[i] = 0;

        while(true)
            {
                o << "T[";

                for(idx_t i = 0; i < (a.rank - 1); ++i)
                    o << ptr[i] << ",";

                o << ptr[N - 1] << "] =";

                o << a.getElem(ptr) << "\n";


                idx_t j;

                for(j = 0; j < N; ++j)
                    {
                        ++ptr[j];

                        if(ptr[j] < a.dimensions[j])
                            break;

                        ptr[j] = 0;
                    }

                if(j == N)
                    break;
            }

        delete [] ptr;

        return o;
    }


    // N = 1 specialization
    /////////////////////////////////////////////////////////////////////////////////////////

  /**
     @brief narray_stream specialization for printing N = 1
   */
    template<typename T>
    struct narray_stream<T, 1>
    {
        narray_stream(void) {}

      //! N = 1 specialized printing operator
        std::ostream &operator()(std::ostream &o, const narray_t<T, 1> &arr);
    };


    template<typename T>
    std::ostream &narray_stream<T, 1>::operator()(std::ostream &o , const narray_t<T, 1> &arr)
    {
        idx_t dim0 = arr.getDim(0);

        o << "[" ;

        for(idx_t i = 0; i < dim0 - 1; ++i)
            o << arr[i] << " , ";

        o << arr[dim0 - 1] << "]\n";

        return o;
    }

    // N = 2 specialization
    /////////////////////////////////////////////////////////////////////////////////////////

  /**
     @brief narray_stream specialization for printing N = 2
   */
    template<typename T>
    struct narray_stream<T, 2>
    {
        narray_stream(void) {}

      //! N = 2 specialized printing operator
        std::ostream &operator()(std::ostream &o, const narray_t<T, 2> &arr);
    };

    template<typename T>
    std::ostream &narray_stream<T, 2>::operator()(std::ostream &o , const narray_t<T, 2> &arr)
    {
        idx_t dim0 = arr.getDim(0);
        idx_t dim1 = arr.getDim(1);

        o << "[";

        for(idx_t row = 0; row < dim0; ++row)
            {
                o << "[";

                for(idx_t col = 0; col < dim1 - 1; ++col)
                    o << arr[row][col] << " , ";

                if(row == dim0 - 1)
                    o << arr[row][dim1 - 1] << "]]\n";
                else
                    o << arr[row][dim1 - 1] << "]\n";
            }

        return o;
    }


    // actual stream operator
    ///////////////////////////////////////////////////////////////////////////////////

  /**
     @brief templated overload on stream operator

     @details calls narray_stream<T,N>::operator() to specialize the printing
     for N = 1 and N = 2 to make things look nice.
   */
    template<typename T, idx_t N>
    std::ostream &operator<<(std::ostream &o, const narray_t<T, N> &a)
    {
        if(!!!a.allocated)
            {
                o << "not allocated, empty tensor" << std::endl;
                return o;
            }

        narray_stream<T, N> ns;
        return ns(o, a);
    }


}
#endif
