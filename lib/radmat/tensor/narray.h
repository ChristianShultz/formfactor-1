#ifndef tensorPRIMITIVE_H_H_GUARD
#define tensorPRIMITIVE_H_H_GUARD





#include <exception>
#include <iostream>
#include <vector>
#include "ensem/ensem.h"
#include "tensor_fwd.h"
#include "radmat/utils/pow2assert.h"

namespace radmat
{

    /**
       @file narray.h
       @brief A collection of templates to generate an n-dimensional array

       @details This headder uses some template recursion to generalize
       an n-dimensional array for simple access using [] operator
       and uses exit codes 101, 102, 103 and 104

    */

    /**
       @brief A recursive template pattern for n-dimensional arrays

       \tparam T the underlying data type
       \tparam N the sub_t 'rank'

       @details a recursive definition for sub arrays slicing with [] operator on each call
       this will allow access as array[d1][d2]...[dn] w/o needing to specialize for
       each possible value of n, basically a generic tensor

       roughly we are mapping a generic n-tensor to linear storage. each application of []
       can be thought of as providing a new view of the tensor as a slice along that index
       the recursive definition of the template operator [] allows us to keep slicing
       untill we hit a single element, the stop point

       this chains together for access, basically in the code you would call something like
       my_array[a][b][c] where my_array is a narray_t<T,3>, the first application of the
       operator [] will leave you with the slice b,c along a and a sub_t<T,2> object.
       my_array<T,3>[a][b][c] -> my_array<T,2>[b][c] -> my_array<T,1>[c] -> tensor element (T^{abc})
     */

    template<typename T, idx_t N>   // type and dimension
    struct sub_t
    {
    public:

        //! access a 'slice' of a sub_t<T,N> returning a sub_t<T,N-1> -- C style index operator
        sub_t < T, N - 1 > operator[](const idx_t index)
        {
            //bound check
            POW2_ASSERT(index < dimensions[0]);
            return sub_t < T, N - 1 > (&elements[index * length[0]], dimensions + 1, length + 1);
        }

        //! C style index operator
        const sub_t < T, N - 1 > operator[](const idx_t index) const
        {
            //bound check
            POW2_ASSERT(index < dimensions[0]);
            return sub_t < T, N - 1 > (&elements[index * length[0]], dimensions + 1, length + 1);
        }

        // declare friends so we can have access to private constructor
        //! grant narray_t access to the private constructor
        friend struct narray_t < T, N + 1 >;
        //! grant sub_t<T,N+1> access to private constructor
        friend struct sub_t < T, N + 1 >;

    private: // data store
        const idx_t *const dimensions;   //! the dimension of this slice
        const idx_t *const length;       //! the number of elements in this slice
        T *const elements;               //! the sub array elements

        //! private constructor so that only other sub_t and narray_t can see it
        sub_t<T, N>(T *_elements, const idx_t *_dimensions, const idx_t *_length)
            : dimensions(_dimensions) , length(_length), elements(_elements)
        {   }

    };



    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////


    /**
       @brief a template partial specialization to specify behavior
       at the recursion stop point, N =1

       \tparam T the underlying data type
       \tparam N the sub_t 'rank'

       @details a recursive definition for sub arrays slicing with [] operator on each call
       this will allow access as array[d1][d2]...[dn] w/o needing to specialize for
       each possible value of n, basically a generic tensor

       roughly we are mapping a generic n-tensor to linear storage. each application of []
       can be thought of as providing a new view of the tensor as a slice along that index
       the recursive definition of the template operator [] allows us to keep slicing
       untill we hit a single element, the stop point

       this chains together for access, basically in the code you would call something like
       my_array[a][b][c] where my_array is a narray_t<T,3>, the first application of the
       operator [] will leave you with the slice b,c along a and a sub_t<T,2> object.
       my_array<T,3>[a][b][c] -> my_array<T,2>[b][c] -> my_array<T,1>[c] -> tensor element (T^{abc})
    */
    template<typename T>
    struct sub_t<T, 1>
    {
        //! C style access to a reference to a number
        T &operator[](const idx_t index)
        {
            //bound check
            POW2_ASSERT(index < dimensions[0]);
            return elements[index];
        }

        //! C syle access
        const T &operator[](const idx_t index) const
        {
            //bound check
            POW2_ASSERT(index < dimensions[0]);
            return elements[index];
        }

        //! grant narray_t<T,1> access to private constructors
        friend struct narray_t<T, 2>;
        //! grant sub_t<T,2> access to private constructor
        friend struct sub_t<T, 2>;

    private:
        const idx_t *const dimensions;   //! the dimension of this slice (the bottom)
        T *const elements;               //! the sub array

        //! private constructor so that only other sub_t and narray_t can see it
        sub_t<T, 1>(T *_elements, const idx_t *_dimensions, const idx_t *_length)
            :  dimensions(_dimensions) , elements(_elements)
        {   }

    };


    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    /**
       @brief the n-dimensional array template

       @details the n-dimensional array type, this will be a data member of the primitive tensor struct.
       this is basically on the top of the sub_t stack and is the container for all elements. Note that sub_t
       never actually took control of anything, it was just a convenient way to define the storage.
       */
    template<typename T, idx_t N>
    struct narray_t
    {
        enum {NN = N};
    public:
        /**
           @brief Initialize a null narray_t of zero length with no elements

           @details Initialize a null narray_t of zero length with no elements
           this is only a data member of a smarter class so no other
           functionality is needed, is allocated by create which will test for successful allocation
        */
        narray_t(void)
            : allocated(false) , rank(NN) , elements(NULL) , nelements(0), dimensions(NULL) , length(NULL)
        {   }

        /**
        @brief Copy Constructor

        @details Create an narray_t from o, allocate new storage and copy elements
         */
        narray_t(const narray_t<T, N> &o)
            : allocated(o.allocated) , rank(o.rank), elements(NULL), nelements(o.nelements) , dimensions(NULL) , length(NULL)
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
           @brief Destructor that cleans up alloacations

           @details we need to explicitly delete here if we allocated
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

        //! return the slice on index -- C style indexing
        sub_t < T, N - 1 > operator[](const idx_t index);

        //! C syle indexing
        const sub_t < T, N - 1 > operator[](const idx_t index) const;

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
        narray_t<T, N> &operator=(const narray_t<T, N> &o);

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

        //! get a 'view' along the elem of idx
        narray_t < T, N - 1 > slice(const idx_t idx, const idx_t elem) const;

        // gcc is stupid, its this way for the friend declaration or the entire
        // function has to be defined in the class, its making a ton of friends..

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
        friend narray_t<T, N> operator-<>(const narray_t<T, N> &plus);

        //! stream operator to make our lives easier higher up
        template<typename U, idx_t M>
        friend std::ostream &operator<<(std::ostream &o, const narray_t<U, M> &arr);

        //! let tensors fiddle about with narray_t pieces as needed
        friend struct Tensor<T, N>;

    private: // data store
        bool allocated;           //! have we allocated storage for this tensor
        idx_t rank;               //! the tensor rank
        T *elements;              //! pointer to array of elements
        idx_t nelements;          //! number of elements
        idx_t *dimensions;        //! pointer to array holding dimensions
        idx_t *length;            //! pointer to this length variable that is convient for mapping the storage
    };



    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////




    // NB this template recursion pattern fails for narray_t<T,1> since the [] operator
    // tries to return a sub_t<T,0> which is ill defined in this framework

    // the fix is to include a full specilization for the narray_t<T,1> struct which is
    // defined in the following include file

} // temporary close of tensor namespace

// include specializations for rank 0 (empty) and rank 1 (works)
#include "narray1.h"

// reopen namespace
namespace radmat
{

    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////

    template<typename T, idx_t N>
    sub_t < T, N - 1 > narray_t<T, N>::operator[](const idx_t index)
    {
        //bound check
        POW2_ASSERT(index < dimensions[0]);
        return sub_t < T, N - 1 > (&elements[index * length[0]], dimensions + 1, length + 1);
    }

    template<typename T, idx_t N>
    const sub_t < T, N - 1 > narray_t<T, N>::operator[](const idx_t index) const
    {
        //bound check
        POW2_ASSERT(index < dimensions[0]);
        return sub_t < T, N - 1 > (&elements[index * length[0]], dimensions + 1, length + 1);
    }

    template<typename T, idx_t N>
    T &narray_t<T, N>::getElem(const std::vector<idx_t> &index)
    {
        idx_t idx(0), count;

        for(count = 0; count < rank; ++count)
            {
                POW2_ASSERT(index[count] < dimensions[count]);  // check that the dimensions make sense
                idx += index[count] * length[count];       // length of the furthest guy is one
            }

        return elements[idx];
    }

    template<typename T, idx_t N>
    const T &narray_t<T, N>::getElem(const std::vector<idx_t> &index) const
    {
        idx_t idx(0), count;

        for(count = 0; count < rank; ++count)
            {
                POW2_ASSERT(index[count] < dimensions[count]);
                idx += index[count] * length[count]; //length of the furthest guy is one
            }

        return elements[idx];
    }

    template<typename T, idx_t N>
    T &narray_t<T, N>::getElem(idx_t const *const *const ptr2ptrs) // assume pointers are of the correct length..
    {
        idx_t idx(0), ct;

        for(ct = 0; ct < rank; ++ct)
            idx += *ptr2ptrs[ct] * length[ct];

        return elements[idx];
    }

    template<typename T, idx_t N>
    const T &narray_t<T, N>::getElem(idx_t const *const *const ptr2ptrs) const // assume pointers are of the correct length..
    {
        idx_t idx(0), ct;

        for(ct = 0; ct < rank; ++ct)
            idx += *ptr2ptrs[ct] * length[ct];

        return elements[idx];
    }

    template<typename T, idx_t N>
    T &narray_t<T, N>::getElem(idx_t const *const ptr)
    {
        idx_t idx(0), ct;

        for(ct = 0; ct < rank; ++ct)
            idx += ptr[ct] * length[ct];

        return elements[idx];
    }

    template<typename T, idx_t N>
    const T &narray_t<T, N>::getElem(idx_t const *const ptr) const
    {
        idx_t idx(0), ct;

        for(ct = 0; ct < rank; ++ct)
            idx += ptr[ct] * length[ct];

        return elements[idx];
    }

    // test allocation to ensure it succeede
    template<typename T, idx_t N>
    bool narray_t<T, N>::create(const idx_t _dimensions[])
    {
        dimensions = new idx_t[rank];
        length = new idx_t[rank];

        nelements = 1;

        for(int i = 0; i < rank; ++i)
            {
                if(_dimensions[i] == 0)
                    {
                        std::cout << "Bad dimensions, allocation failing" << std::endl;
                        return false;
                    }

                nelements *= _dimensions[i];
                dimensions[i] = _dimensions[i];

                length[i] = 1;

                for(int j = N - 1; j > i; --j)
                    length[i] *= _dimensions[j];  // this has to be _dimension since dimension
            }                                 // dimension has yet to be filled fully.. that

        // was an annoying typo to track down, sigh
        elements = new T[nelements];

        // bad allocation
        if(!elements)
            {
                std::cout << "Could not allocate storage" << std::endl;
                return false;
            }

	  for(idx_t elem = 0; elem < nelements; ++elem)
            elements[elem] = 0.;

        allocated = true;
        return true;
    }

    // copy assignment
    template<typename T, idx_t N>
    narray_t<T, N> &narray_t<T, N>::operator=(const narray_t<T, N> &o)
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

    template<typename T, idx_t N>
    narray_t < T, N - 1 > narray_t<T, N>::slice(const idx_t _idx, const idx_t elem) const
    {
        // allocate and set up indicies for looping
        idx_t *cycle, **slice_, **full, el(elem);
        bool minus = false;
        std::vector<idx_t> dim(rank - 1, 0);

        cycle = new idx_t[N - 1];
        slice_ = new idx_t*[N - 1];
        full = new idx_t*[N];

        for(idx_t i = 0; i < rank; ++i)
            {
                if(i < rank - 1)
                    {
                        cycle[i] = 0;
                        slice_[i] = &cycle[i];
                    }

                if(i != _idx)
                    if(!minus)
                        {
                            full[i] = &cycle[i];
                            dim[i] = dimensions[i];
                        }
                    else
                        {
                            full[i] = &cycle[i - 1];
                            dim[i - 1] = dimensions[i];
                        }
                else
                    {
                        minus = true;
                        full[i] = &el;
                    }
            }

        narray_t < T, N - 1 > ret;

        if(!ret.create(&dim[0]))
            {
                std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
                std::cout << "Allocation error during tensor slicing" << std::endl;
                exit(103);
            }

        const idx_t sz = dim.size();

        while(true)
            {
                ret.getElem(slice_) = this->getElem(full);
                idx_t j;

                for(j = 0; j < sz; ++j)
                    {
                        ++cycle[j];

                        if(cycle[j] < dim[j])
                            break;

                        cycle[j] = 0;
                    }

                if(j == N - 1)
                    break;
            }

        delete [] cycle;
        delete [] slice_;
        delete [] full;

        return ret;
    }


    template<typename T, typename U, idx_t N, idx_t M >
    narray_t < typename ENSEM::Promote<T, U>::Type_t, N + M - 2 > contract(const narray_t<T, N> &A,
            const narray_t<U, M> &B,
            const idx_t A_idx,
            const idx_t B_idx)
    {
        // check the dimensions
        POW2_ASSERT(A.dimensions[A_idx] == B.dimensions[B_idx]);

        const idx_t c_dim = A.dimensions[A_idx];    // dimension of the index we are contracting over
        idx_t contract_idx;                         // the contraction index

        // allocate some storage
        idx_t **a, **b, *c;
        a = new idx_t*[N];                          // specifies the index in the A tensor for getElem
        b = new idx_t*[M];                          // same for B
        c = new idx_t[N + M - 2];                   // same for the return tensor

        // set up the pointer arrays for A
        bool minus = false;                         // c will specify the index in the new array

        for(idx_t i = 0; i < A.rank; ++i)           // since it is a contraction we can use the
            if(i != A_idx)                            // bits of c to tell us what elements we needed
                if(!minus)                              // from a and b, using pointers to the c array
                    a[i] = &c[i];                         // means that we don't need to worry about
                else                                    // keeping the indicies together.
                    a[i] = &c[i - 1];
            else                                      // the set of lines following 'allocate some storage'
                {
                    // declare the storage for the c array which is the
                    a[i] = &contract_idx;                 // index in the new tensor and then tie the pointers
                    minus = true;                         // for a and b to their corresponding index in c
                }

        // set up the pointer arrays for B
        minus = false;

        for(idx_t i = 0; i < B.rank; ++i)
            if(i != B_idx)
                if(!minus)
                    b[i] = &c[i + A.rank - 1];
                else
                    b[i] = &c[i + A.rank - 2];
            else
                {
                    b[i] = &contract_idx;
                    minus = true;
                }

        // set up the array to be returned
        narray_t < typename ENSEM::Promote<T, U>::Type_t, N + M - 2 > ret;

        // get the dimensions
        std::vector<idx_t> dim;

        for(idx_t i = 0; i < A.rank; ++i)
            if(i != A_idx)
                dim.push_back(A.dimensions[i]);

        for(idx_t i = 0; i < B.rank; ++i)
            if(i != B_idx)
                dim.push_back(B.dimensions[i]);

        // allocate the storage
        if(!ret.create(&dim[0]))
            {
                std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n"
                          << "Allocation error during contraction" << std::endl;
                exit(102);
            }

        // Now fill the return array

        // zero the index explicitly
        for(idx_t i = 0; i < A.rank + B.rank - 2; ++i)
            c[i] = 0;

        const idx_t sz = dim.size();

        // set up an n-dimensional for loop to loop over an fill in all the slots
        while(true)
            {
                // set the element to zero
                ret.getElem(c) = 0;

                // a,b contain pointers to contract_idx so looping over this is like saying
                // T^{abc} = \sum_e T^{ae}*T^{ebc}
                for(contract_idx = 0; contract_idx < c_dim; ++contract_idx)
                    ret.getElem(c) += A.getElem(a) * B.getElem(b);

                // update idx
                idx_t j;

                for(j = 0; j < sz; ++j)        // the leftmost index cycles fastest
                    {
                        // perhaps it should be changed to be
                        ++c[j];                            // the rightmost based on how the arrays

                        if(c[j] < dim[j])                  // are stored internally?
                            break;

                        c[j] = 0;
                    }

                // break while loop
                if(j == dim.size())
                    break;
            }

        // deallocate
        delete [] c;     // this was an array of idx_t
        delete [] a;     // this was an array of pointers to idx_t type
        delete [] b;     // same as a

        return ret;
    }


    template<typename T, typename U, idx_t M>
    narray_t<typename ENSEM::Promote<T, U>::Type_t, M> contract_M(const narray_t<T, M> &A, const narray_t<U, 2> &B, const idx_t A_idx)
    {
        POW2_ASSERT(A.dimensions[A_idx] == B.dimensions[0]);
        POW2_ASSERT(B.dimensions[0] == B.dimensions[1]);
        POW2_ASSERT(A_idx < A.rank);

        idx_t **a, **b, *cycle, contract;
        std::vector<idx_t> dim;


        a = new idx_t*[A.rank];
        b = new idx_t*[2];
        cycle = new idx_t[A.rank];

        for(idx_t i = 0; i < A.rank; ++i)
            {
                dim.push_back(A.dimensions[i]);

                if(i != A_idx)
                    a[i] = &cycle[i];
                else
                    a[i] = &contract;
            }

        b[0] = &contract;
        b[1] = &cycle[A_idx];

        narray_t<typename ENSEM::Promote<T, U>::Type_t, M> ret;

        if(!ret.create(&dim[0]))
            {
                std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n"
                          << "Allocation error during contraction" << std::endl;
                exit(104);
            }

        const idx_t cdim = A.dimensions[A_idx];
        const idx_t sz = dim.size();

        while(true)
            {
                ret.getElem(cycle) = 0.;

                for(contract = 0; contract < cdim; ++contract)
                    ret.getElem(cycle) += A.getElem(a) * B.getElem(b);

                idx_t j;

                for(j = 0; j < sz; ++j)
                    {
                        ++cycle[j];

                        if(cycle[j] < dim[j])
                            break;

                        cycle[j] = 0;
                    }

                if(j == sz)
                    break;
            }

        delete [] cycle;
        delete [] a;
        delete [] b;

        return ret;
    }

    template<typename T, typename U, idx_t M>
    narray_t<typename ENSEM::Promote<T, U>::Type_t, M> operator+(const narray_t<T, M> &lhs, const narray_t<U, M> &rhs)
    {
        std::vector<idx_t> dim(M, 0);

        for(idx_t idx = 0; idx < M; ++idx)
            {
                dim[idx] = lhs.dimensions[idx];
                POW2_ASSERT(lhs.dimensions[idx] == rhs.dimensions[idx]);
            }

        narray_t<typename ENSEM::Promote<T, U>::Type_t, M> ret;
        ret.create(&dim[0]);
        const idx_t nelem = lhs.nelements;

        for(idx_t idx = 0; idx < nelem; ++idx)
            ret.elements[idx] = lhs.elements[idx] + rhs.elements[idx];

        return ret;
    }


    template<typename T, typename U, idx_t M>
    narray_t<typename ENSEM::Promote<T, U>::Type_t, M> operator*(const narray_t<T, M> &lhs,  const U &rhs)
    {
        std::vector<idx_t> dim(M, 0);

        for(idx_t idx = 0; idx < M; ++idx)
            dim[idx] = lhs.dimensions[idx];

        narray_t<typename ENSEM::Promote<T, U>::Type_t, M> ret;
        ret.create(&dim[0]);
        const idx_t nelem = lhs.nelements;

        for(idx_t idx = 0; idx < nelem; ++idx)
            ret.elements[idx] = lhs.elements[idx] * rhs;

        return ret;
    }

    template<typename T, typename U, idx_t M>
    narray_t<typename ENSEM::Promote<T, U>::Type_t, M> operator*(const U &rhs, const narray_t<T, M> &lhs)
    {
        return lhs * rhs;
    }

    template<typename T, typename U, idx_t M>
    narray_t<typename ENSEM::Promote<T, U>::Type_t, M> operator/(const narray_t<T, M> &lhs,  const U &rhs)
    {
        return lhs * (T(1.) / rhs);
    }

    template<typename T, typename U, idx_t M>
    narray_t<typename ENSEM::Promote<T, U>::Type_t, M> operator-(const narray_t<T, M> &lhs,  const narray_t<U, M> &rhs)
    {
        narray_t<U, M> dum = -rhs;
        return lhs + dum;
    }


    template<typename T, idx_t N>
    narray_t<T, N> operator-(const narray_t<T, N> &plus)
    {
        narray_t<T, N> minus;

        minus.allocated = plus.allocated;
        minus.rank = plus.rank;
        minus.nelements = plus.nelements;

        minus.dimensions = new idx_t[N];
        minus.length = new idx_t[N];
        minus.elements = new T[minus.nelements];

        memcpy(minus.dimensions, plus.dimensions, N * sizeof(idx_t));
        memcpy(minus.length, plus.length, N * sizeof(idx_t));

        for(idx_t i = 0; i < minus.nelements; ++i)
            minus.elements[i] = -plus.elements[i];

        return minus;
    }


} // close tensor namespace



#endif
