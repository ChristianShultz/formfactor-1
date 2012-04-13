#ifndef NARRAY1_H_H_GUARD
#define NARRAY1_H_H_GUARD

#include <exception>
#include <iostream>
#include <vector>
#include "ensem/ensem.h"
#include "tensor_fwd.h"
#include "utils/pow2assert.h"


/**
   @file narray1.h
   @brief Specializations for narray<T,1> and narray<T,0>

   @details this template recursion pattern fails for narray_t<T,1> since the [] operator
   tries to return a sub_t<T,0> which is ill defined in this framework

   the fix is to include a full specilization for the narray_t<T,1> struct which is
   defined in this include file
 */

namespace tensor
{
    /**
       @brief an empty specialization to prevent improper use
     */

    template<typename T>
    struct narray_t<T, 0>
    {

    };

    /**
       @brief Specialization for narray<T,1> struct

       @details this template recursion pattern fails for narray_t<T,1> since the [] operator
       tries to return a sub_t<T,0> which is ill defined in this framework

       the fix is to include a full specilization for the narray_t<T,1> struct which is
       defined in this include file
    */
    template<typename T>
    struct narray_t<T, 1>
    {

        enum {NN = 1};
    public:
        /** @brief initialize a null narray_t of zero length with no elements
          @details this is only a data member of a smarter class so no other
           functionality is needed, is allocated by bool create or
           by void alloc, create will test for successful allocation
        */
        narray_t(void)
            : allocated(false) , rank(NN) , elements(NULL) , nelements(0), dimensions(NULL) , length(NULL)
        {   }

        /**
        @brief Copy Constructor

        @details Create an narray_t from o, allocate new storage and copy elements
         */
        narray_t(const narray_t<T, 1> &o)
            : allocated(o.allocated) , 
	      rank(o.rank), 
	      elements(NULL), 
	      nelements(o.nelements) , 
	      dimensions(NULL) , 
	      length(NULL)
        {
            if(allocated)
                {
                    dimensions = new idx_t[rank];
                    length = new idx_t[rank];
                    elements = new T[nelements];

                    // this should be faster than a for loop
                    memcpy(dimensions, o.dimensions, rank * sizeof(idx_t));
                    memcpy(length, o.length, rank * sizeof(idx_t));
		    // memcpy(elements, o.elements, nelements * sizeof(T));

		    // what if T is something strange -- use copy assignment
		    for(idx_t i = 0; i < nelements; ++i)
		      elements[i] = o.elements[i];
                }
        }

        /**
        @brief Copy Constructor

        @details Create an narray_t from o, allocate new storage and copy elements
         */
        ~narray_t(void)
        {
            if(allocated)
                {
                    delete [] elements;
                    delete [] dimensions;
                    delete [] length;
                }
        }

        //! return a ref to a number -- C style indexing
        T &operator[](const idx_t index);

        //! C style indexing
        const T &operator[](const idx_t index) const;

        /** @brief a 'hack' to allow access to an element by a vector of indicies
        @details a 'hack' to allow for access by a vector of indicies which is nice for looping in
        an arbitrary number of dimensions as we surely must at some point
         */
        T &getElem(const std::vector<idx_t> &index);

        /** @brief a 'hack' to allow access to an element by a vector of indicies
        @details a 'hack' to allow for access by a vector of indicies which is nice for looping in
        an arbitrary number of dimensions as we surely must at some point
         */
        const T &getElem(const std::vector<idx_t> &index) const;


        //! contraction helper, no safety checks on inputs
        T &getElem(idx_t const *const *const ptr2ptrs);

        //! contraction helper, no safety checks on inputs
        const T &getElem(idx_t const *const *const ptr2ptrs) const;

        //! contraction helper, no safety checks on inputs
        T &getElem(idx_t const *const ptr);

        //! contraction helper, no safety checks on inputs
        const T &getElem(idx_t const *const ptr) const;

        //! test allocation to ensure it succeede
        bool create(const idx_t _dimensions[]);

        //! copy assignment
        narray_t<T, 1> &operator=(const narray_t<T, 1> &o);

        idx_t getRank(void) const
        {
            return rank;
        }

        //! return dimension on index i
        idx_t getDim(const idx_t i) const
        {
            POW2_ASSERT(i < rank);
            return dimensions[i];
        }

        //! doesn't actually work, asserts false
        sub_t <T, 0> slice(const idx_t idx, const idx_t elem) const
        {
            POW2_ASSERT(false);  // what were you doing calling this function?
        }

        /**
        @brief this contraction takes \f$ T^{a..b..c} \otimes T^{p..b..q} \rightarrow T^{a..c,p..q} \f$

        @details contract two arrays along the specified index to return this nasty thing
         this is the main workhorse so if someone smarter than me wants to figure out
         how to write this more efficiently I would be much obliged.

         this contraction takes \f$ T^{a..b..c} \otimes T^{p..b..q} \rightarrow T^{a..c,p..q} \f$
        */
        template<typename U, typename UU, idx_t M, idx_t MM>
        friend narray_t < typename ENSEM::Promote<U, UU>::Type_t, M + MM - 2 > contract(const narray_t<U, M> &A,
                const narray_t<UU, MM> &B,
                const idx_t A_idx,
                const idx_t B_idx);

        /**
           @brief this contraction is for metrics \f$ T^{a..b..c}  g_{bq} \rightarrow T^{a.. \:\; ..c}_{\quad q} \f$

           @details this contraction is for metric and leaves the index in place.
           it takes \f$ T^{a..b..c}  g_{bq} \rightarrow T^{a.. \:\; ..c}_{\quad q} \f$ , it exists
           to raise/lower w/o disturbing the order of the indicies and is a special
           case of the general contractions
        */
        template<typename U, typename UU, idx_t M>
        friend narray_t <typename ENSEM::Promote<U, UU>::Type_t, M> contract_M(const narray_t<U, M> &A,
                const narray_t<UU, 2> &B,
                const idx_t A_idx);

        //! algebraic operator to make our lives easier higher up
        template<typename U, typename UU, idx_t M>
        friend narray_t <typename ENSEM::Promote<U, UU>::Type_t, M> operator+(const narray_t<U, M> &lhs, const narray_t<UU, M> &rhs);

        //! algebraic operator to make our lives easier higher up
        template<typename U, typename UU, idx_t M>
        friend narray_t <typename ENSEM::Promote<U, UU>::Type_t, M > operator*(const narray_t<U, M> &lhs, const UU &rhs);

        //! algebraic operator to make our lives easier higher up
        template<typename U, typename UU, idx_t M>
        friend narray_t <typename ENSEM::Promote<U, UU>::Type_t, M> operator/(const narray_t<U, M> &lhs, const UU &rhs);

        //! algebraic operator to make our lives easier higher up
        template<typename U, typename UU, idx_t M>
        friend narray_t<typename ENSEM::Promote<U, UU>::Type_t, M> operator-(const narray_t<U, M> &lhs, const narray_t<UU, M> &rhs);

        //! algebraic operator to make our lives easier higher up
        friend narray_t<T, 1> operator-<>(const narray_t<T, 1> &plus);

        //! stream operator to make our lives easier higher up
        template<typename U, idx_t M>
        friend std::ostream &operator<<(std::ostream &o, const narray_t<U, M> &arr);

        //! let tensors fiddle about with narray_t pieces as needed
        friend struct Tensor<T, 1>;

    private: // data store
        bool allocated;           //! have we allocated storage for this tensor
        idx_t rank;               //! the tensor rank
        T *elements;              //! pointer to array of elements
        idx_t nelements;          //! number of elements
        idx_t *dimensions;        //! pointer to array holding dimensions
        idx_t *length;            //! pointer to this length variable that is convient for mapping the storage
    };


    template<typename T>
    T &narray_t<T, 1>::operator[](const idx_t index)
    {
        //bound check
        POW2_ASSERT(index < dimensions[0]);
        return elements[index];
    }

    template<typename T>
    const T &narray_t<T, 1>::operator[](const idx_t index) const
    {
        //bound check
        POW2_ASSERT(index < dimensions[0]);
        return elements[index];
    }

    template<typename T>
    T &narray_t<T, 1>::getElem(const std::vector<idx_t> &index)
    {
        return (*this)[index[0]];
    }

    template<typename T>
    const T &narray_t<T, 1>::getElem(const std::vector<idx_t> &index) const
    {
        return (*this)[index[0]];
    }

    template<typename T>
    T &narray_t<T, 1>::getElem(idx_t const *const *const ptr2ptr)
    {
        return (*this)[*ptr2ptr[0]];
    }

    template<typename T>
    const T &narray_t<T, 1>::getElem(idx_t const *const *const ptr2ptr) const
    {
        return (*this)[*ptr2ptr[0]];
    }

    template<typename T>
    T &narray_t<T, 1>::getElem(idx_t  const *const ptr2ptr)
    {
        return (*this)[ptr2ptr[0]];
    }

    template<typename T>
    const T &narray_t<T, 1>::getElem(idx_t const *const ptr2ptr) const
    {
        return (*this)[ptr2ptr[0]];
    }

    template<typename T>
    bool narray_t<T, 1>::create(const idx_t _dimensions[])
    {
        dimensions = new idx_t[1];
        length = new idx_t[1];
        nelements = _dimensions[0];
        length[0] = _dimensions[0];
        dimensions[0] = _dimensions[0];
        elements = new T[nelements];

        POW2_ASSERT(dimensions);
        POW2_ASSERT(length);
        POW2_ASSERT(elements);

        for(idx_t elem = 0; elem < nelements; ++elem)
	  elements[elem] = 0.;

        allocated = true;

        return true;
    }


    // copy assignment
    template<typename T>
    narray_t<T, 1> &narray_t<T, 1>::operator=(const narray_t<T, 1> &o)
    {
        if(this != &o)
            {
                rank = o.rank;
                nelements = o.nelements;

                if(allocated)
                    {
                        delete [] elements;
                        delete [] dimensions;
                        delete [] length;
                    }
                else
                    allocated = true;

                dimensions = new idx_t[rank];
                length = new idx_t[rank];
                elements = new T[nelements];

                // this should be faster than a for loop
                memcpy(dimensions, o.dimensions, rank * sizeof(idx_t));
                memcpy(length, o.length, rank * sizeof(idx_t));
		    // memcpy(elements, o.elements, nelements * sizeof(T));

		    // what if T is something strange -- use copy assignment
		    for(idx_t i = 0; i < nelements; ++i)
		      elements[i] = o.elements[i];
            }

        return *this;
    }

}
#endif
