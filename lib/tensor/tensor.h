#ifndef tensor_H_H_GUARD
#define tensor_H_H_GUARD

#include <iostream>
#include "tensor_fwd.h"
#include "narray.h"
#include "narray_stream.h"
#include "ensem/ensem.h"
#include "utils/pow2assert.h"

// this file uses exit codes 201

namespace tensor
{
  /**
     @file tensor.h
     @brief the template for rank N tensors
     
     @details this file uses exit codes 201
   */

  /**
     @brief A pure virtual base class to derive from

     @details The factory class for polarisation vectors will store pointers to 
     this base class so that we only need one factory for all ranks
   */
    struct TensorImplBase
    {
        virtual TensorImplBase *clone(void) const
        {
            return new TensorImplBase(*this);
        }
    };

  /** 
      @brief some meta programming magic to specify behavior for rank 1 tensors.. sigh
   */
    template<typename T, idx_t N>
    struct subRet
    {
        typedef sub_t < T, N - 1 > Type_t;
    };

  /** 
      @brief some meta programming magic to specify behavior for rank 1 tensors.. sigh
   */
    template<typename T>
    struct subRet<T, 1>
    {
        typedef T &Type_t;
    };

  /** 
      @brief some meta programming magic to specify behavior for rank 1 tensors.. sigh
   */
    template<typename T, idx_t N>
    struct subRetc
    {
        typedef const sub_t < T, N - 1 > Type_t;
    };

  /** 
      @brief some meta programming magic to specify behavior for rank 1 tensors.. sigh
   */
    template<typename T>
    struct subRetc<T, 1>
    {
        typedef const T &Type_t;
    };


  /**
     @brief The rank N tensor template

     @details Holds a generic narray_t and a vector of bools to keep track of the position
     of indicies for contraction, derives from TensorImplBase to allow for polymorphism 
   */
    template<typename T, idx_t N>
    struct Tensor : public TensorImplBase
    {
    public: // constructor destructor assignment
      //! creates an empty rank N tensor with all indicies raised -- set to true
        Tensor(void);

      //! creates a rank N tensor from _data with all indicies raised  -- set to true
        Tensor(const narray_t<T, N> &_data);

      //! copy constructor
        Tensor(const Tensor<T, N> &o);

      //! assignment operator
        Tensor<T, N> &operator=(const Tensor<T, N> &o);

      //! do nothing destructor
        ~Tensor<T, N>(void) {}

      //! create the tensor with the appropriate dimensions
        bool create(const std::vector<idx_t> &dimensions);

      //! create the tensor with the appropriate dimensions
        bool create(const idx_t dimensions[]);

      /**
	 @brief clone a tensor returning a new instance

	 @details this allocates a new object, it is the responsibility of the
	 user to ensure that it is properly deallocated somewhere else in the 
	 program
       */
        Tensor<T, N> *clone(void)
        {
            return new Tensor<T, N>(*this);
        }

      //! C style access operator
        typename subRet < T, N >::Type_t operator[](const idx_t index);

      //! C style access operator
        typename subRetc< T, N >::Type_t operator[](const idx_t index) const;

      //! access a element using a vector of indicies
        T &getElem(const std::vector<idx_t> &index);

      //! access a element using a vector of indicies
        const T &getElem(const std::vector<idx_t> &index) const;

      //! access an element, doesn't check length of pointer/dimensions
        T &getElem(const idx_t *const ptr);       

      //! access an element, doesn't check length of pointer/dimensions                     
        const T &getElem(const idx_t *const ptr) const;   

      //! access an element, doesn't check length of pointer/dimensions                 
        T &getElem(const idx_t *const *const ptr2ptrs);            

      //! access an element, doesn't check length of pointer/dimensions         
        const T &getElem(const idx_t *const *const ptr2ptrs) const;          


        // index

      //! check to see if an index is raised
        bool isRaised(const idx_t idx) const;    

      //! lower index w/o application of a metric                         
        void lowerNoMetric(const idx_t index);

      //! raise index w/o application of a metric                              
        void raiseNoMetric(const idx_t index);

      //! raise index by application of metric
        void raiseApplyMetric(const Tensor<T, 2> &metric, const idx_t idx); 

      //! lower index by application of a metric
        void lowerApplyMetric(const Tensor<T, 2> &metric, const idx_t idx);

      //! bulk set the up (down) position of indicies with a vector
        void setIndicies(const std::vector<bool> &idxs);                   

      //! flip all the indicies without application of a metric
        void flipAllIndicies(void);

      //! flip all the indicies with application of a metric
        void flipAllIndicies(const Tensor<T, 2> &metric);

      //! return the rank of the tensor
        idx_t getRank(void) const
        {
            return tensor.getRank();
        }

      //! return the dimension of an index
        idx_t getDim(const idx_t i) const
        {
            return tensor.getDim(i);
        }

        //! algebraic operator
        Tensor<T, N> &operator+=(const Tensor<T, N> &o);

        //! algebraic operator
        Tensor<T, N> &operator+=(const T some_const);

        //! algebraic operator
        Tensor<T, N> &operator-=(const Tensor<T, N> &o);

        //! algebraic operator
        Tensor<T, N> &operator-=(const T some_const);

        //! algebraic operator
        Tensor<T, N> &operator*=(const T some_const);

        //! algebraic operator
        Tensor<T, N> &operator/=(const T some_const);

        //! slice a tensor providing a 'view' along the index at the element
        Tensor < T, N - 1 > slice(const idx_t idx, const idx_t elem) const;

        //! help with slicing
        friend struct Tensor < T, N + 1 >;

        //! unary minus
        friend Tensor<T, N> operator-<>(const Tensor<T, N> &plus);


      /**
	 @brief raise (lower) an index by application of a metric

	 @ details ignores indicies on the metric and raises (lowers)
	 the index of the idx of the return tensor after contracting on the primitive narray_t types
       */
        template<typename U, typename UU, idx_t M>
        friend Tensor<typename ENSEM::Promote<U, UU>::Type_t, M> applyMetric(const Tensor<U, M> &_tensor,
                const Tensor<UU, 2> &metric,
                const idx_t idx);

        //! contract without a metric, raised/lowered checking, can not be used in place of applyMetric
        template<typename U, typename UU, idx_t M, idx_t MM>
        friend Tensor < typename ENSEM::Promote<U, UU>::Type_t, M + MM - 2 > contract(const Tensor<U, M> &A,
                const Tensor<UU, MM> &B,
                const idx_t A_idx,
                const idx_t B_idx);

        //! contract with a metric, apply the metric if needed
        template<typename U, typename UU, typename UUU, idx_t M, idx_t MM>
        friend Tensor < typename ENSEM::Promote<typename ENSEM::Promote<U, UU>::Type_t, UUU>::Type_t, M + MM - 2 > contract(const Tensor<U, M> &A,
                const Tensor<UU, MM> &B,
                const Tensor<UUU, 2> &metric,
                const idx_t A_idx,
                const idx_t B_idx);

        //! basic asymmetric algebraic operation
        template<typename U, typename UU, idx_t M>
        friend Tensor<typename ENSEM::Promote<U, UU>::Type_t, M> operator*(const Tensor<U, M> &lhs, const UU rhs);

        //! basic asymmetric algebraic operation
        template<typename U, typename UU, idx_t M>
        friend Tensor<typename ENSEM::Promote<U, UU>::Type_t, M> operator/(const Tensor<U, M> &lhs, const UU rhs);

      /**
	 @brief basic asymmetric algebraic operation
	 @details indicies of lhs and rhs must be in same positions for this operation to make sense
      */
        template<typename U, typename UU, idx_t M>
        friend Tensor<typename ENSEM::Promote<U, UU>::Type_t, M> operator+(const Tensor<U, M> &lhs, const Tensor<UU, M> &rhs);

      /**
	 @brief basic asymmetric algebraic operation
	 @details indicies of lhs and rhs must be in same positions for this operation to make sense
      */
        template<typename U, typename UU, idx_t M>
        friend Tensor<typename ENSEM::Promote<U, UU>::Type_t, M> operator-(const Tensor<U, M> &lhs, const Tensor<UU, M> &rhs);

        //! stream operator
        template<typename U, idx_t M>
        friend std::ostream &operator<<(std::ostream &o, const Tensor<U, M> &t);

    private: // data store
      narray_t<T, N> tensor;      //! the elements of the tensor 
      std::vector<bool> raised;   //! the position of the indicies, true is raised
    };


    // Implementation
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////

    template<typename T, idx_t N>
    Tensor<T, N>::Tensor(void)
    {
        raised.resize(N, true);
    }

    template<typename T, idx_t N>
    Tensor<T, N>::Tensor(const narray_t<T, N> &_data)
        : tensor(_data)
    {
        raised.resize(N, true);
    }

    template<typename T, idx_t N>
    Tensor<T, N>::Tensor(const Tensor<T, N> &o)
    {
        tensor = o.tensor;
        raised = o.raised;
    }

    template<typename T, idx_t N>
    Tensor<T, N> &Tensor<T, N>::operator=(const Tensor<T, N> &o)
    {
        if(this != &o)
            {
                tensor = o.tensor;
                raised = o.raised;
            }

        return *this;
    }

    template<typename T, idx_t N>
    bool Tensor<T, N>::create(const std::vector<idx_t> &dimensions)
    {
        raised.resize(N, true);
        return tensor.create(&dimensions[0]);
    }

    template<typename T, idx_t N>
    bool Tensor<T, N>::create(const idx_t dimensions[])
    {
        raised.resize(N, true);
        return tensor.create(dimensions);
    }

    template<typename T, idx_t N>
    typename subRet<T, N>::Type_t Tensor<T, N>::operator[](const idx_t index)
    {
        return tensor[index];
    }

    template<typename T, idx_t N>
    typename subRetc<T, N>::Type_t Tensor<T, N>::operator[](const idx_t index) const
    {
        return tensor[index];
    }

    template<typename T, idx_t N>
    T &Tensor<T, N>::getElem(const std::vector<idx_t> &idx)
    {
        return tensor.getElem(idx);
    }

    template<typename T, idx_t N>
    const T &Tensor<T, N>::getElem(const std::vector<idx_t> &idx) const
    {
        return tensor.getElem(idx);
    }

    template<typename T, idx_t N>
    T &Tensor<T, N>::getElem(const idx_t *const ptr)
    {
        return tensor.getElem(ptr);
    }

    template<typename T, idx_t N>
    const T &Tensor<T, N>::getElem(const idx_t *const ptr) const
    {
        return tensor.getElem(ptr);
    }

    template<typename T, idx_t N>
    T &Tensor<T, N>::getElem(const idx_t *const *const ptr2ptrs)
    {
        return tensor.getElem(ptr2ptrs);
    }

    template<typename T, idx_t N>
    const T &Tensor<T, N>::getElem(const idx_t *const *const ptr2ptrs) const
    {
        return tensor.getElem(ptr2ptrs);
    }

    template<typename T, idx_t N>
    bool Tensor<T, N>::isRaised(const idx_t idx) const
    {
        POW2_ASSERT(idx < raised.size());
        return raised[idx];
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::lowerNoMetric(const idx_t idx)
    {
        POW2_ASSERT(idx < raised.size());

        if(!raised[idx])
            {
                std::cout << "Warning: " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n"
                          << "lowering an already lowered index" << std::endl;
            }

        raised[idx] = false;
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::raiseNoMetric(const idx_t idx)
    {
        POW2_ASSERT(idx < raised.size());

        if(raised[idx])
            {
                std::cout << "Warning: " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n"
                          << "raising an already raised index" << std::endl;
            }

        raised[idx] = true;
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::raiseApplyMetric(const Tensor<T, 2> &metric, const idx_t idx)
    {
        if(raised[idx])
            {
                std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
                std::cout << "Raising an index by application of metric on an already raised index, exiting" << std::endl;
                exit(201);
            }

        tensor = contract_M(tensor, metric, idx);
        raised[idx] = true;
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::lowerApplyMetric(const Tensor<T, 2> &metric, const idx_t idx)
    {
        if(!raised[idx])
            {
                std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
                std::cout << "Lowering an index by application of metric on an already lowered index, exiting" << std::endl;
                exit(202);
            }

        tensor = contract_M(tensor, metric, idx);
        raised[idx] = false;
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::setIndicies(const std::vector<bool> &idxs)
    {
        POW2_ASSERT(idxs.size() == raised.size());
        raised = idxs;
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::flipAllIndicies(void)
    {
        std::vector<bool>::iterator it;

        for(it = raised.begin(); it != raised.end(); ++it)
            *it = !*it;
    }

    template<typename T, idx_t N>
    void Tensor<T, N>::flipAllIndicies(const Tensor<T, 2> &metric)
    {
        for(idx_t i = 0; i < raised.size(); ++i)
            if(raised[i])
                lowerApplyMetric(metric, i);
            else
                raiseApplyMetric(metric, i);
    }


    template<typename T, idx_t N>
    Tensor<T, N>  &Tensor<T, N>::operator+=(const Tensor<T, N> &o)
    {
        // check index
        for(idx_t i = 0; i < N; ++i)
            POW2_ASSERT(raised[i] == o.raised[i]);

        tensor = (tensor + o.tensor);

        return *this;
    }

    template<typename T, idx_t N>
    Tensor<T, N> &Tensor<T, N>::operator+=(const T some_const)
    {
        for(idx_t idx = 0; idx < tensor.nelements; ++idx)
            tensor.elements[idx] += some_const;

        return *this;
    }

    template<typename T, idx_t N>
    Tensor<T, N> &Tensor<T, N>::operator-=(const Tensor<T, N> &o)
    {
        return (*this) += (-o);
    }

    template<typename T, idx_t N>
    Tensor<T, N> &Tensor<T, N>::operator-=(const T some_const)
    {
        return (*this) += (-some_const);
    }

    template<typename T, idx_t N>
    Tensor<T, N> &Tensor<T, N>::operator*=(const T some_const)
    {
        for(idx_t i = 0; i < tensor.nelements; ++i)
            tensor.elements[i] *= some_const;

        return *this;
    }

    template<typename T, idx_t N>
    Tensor<T, N> &Tensor<T, N>::operator/=(const T some_const)
    {
        for(idx_t i = 0; i < tensor.nelements; ++i)
            tensor.elements[i] /= some_const;

        return *this;
    }

    template<typename T, idx_t N>
    Tensor < T, N - 1 > Tensor<T, N>::slice(const idx_t idx, const idx_t elem) const
    {
        narray_t < T, N - 1 > tmp = tensor.slice(idx, elem);
        Tensor < T, N - 1 > ret(tmp);
        std::vector<bool> isRaised(raised);

        isRaised.erase(isRaised.begin() + idx);
        ret.raised = isRaised;

        return ret;
    }


    // Friends
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    // unary minus
    template<typename T, idx_t N>
    Tensor<T, N> operator-(const Tensor<T, N> &plus)
    {
        Tensor<T, N> minus(plus);

        minus.tensor = (-plus.tensor);

        return minus;
    }


    // apply a metric to raise (lower) an index
    template<typename T, typename U, idx_t M>
    Tensor<typename ENSEM::Promote<T, U>::Type_t, M> applyMetric(const Tensor<T, M> &_tensor,
            const Tensor<U, 2> &metric,
            const idx_t idx)
    {

        if(_tensor.isRaised(idx) == metric.isRaised(0))
            {
                std::cout << "Warning: " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
                std::cout << "metric index and tensor index are both raised (lowered)," <<
                          " proceeding as if they weren't" << std::endl;
            }

        narray_t<T, M> tmp = contract_M(_tensor.tensor, metric.tensor, idx);

        Tensor<typename ENSEM::Promote<T, U>::Type_t, M> ret = _tensor;
        ret.tensor = tmp;

        if(ret.isRaised(idx))
            ret.raised[idx] = false;
        else
            ret.raised[idx] = true;

        return ret;
    }

    template<typename T, typename U, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<T, U>::Type_t , N + M - 2 > contract(const Tensor<T, N> &A,
            const Tensor<U, M> &B,
            const idx_t A_idx,
            const idx_t B_idx)
    {
        // contract on primitive types
        POW2_ASSERT(A.isRaised(A_idx) != B.isRaised(B_idx));
        narray_t < typename ENSEM::Promote<T, U>::Type_t, M + N - 2 > ret_tens;
        ret_tens = contract(A.tensor, B.tensor, A_idx, B_idx);
        Tensor < typename ENSEM::Promote<T, U>::Type_t, M + N - 2 > ret;
        ret.tensor = ret_tens;
        std::vector<bool> isRaised;

        // determine the positions of the indicies
        // !!! NB: this relies on how the narray_t contractions work to get it
        // right so be careful if things change
        for(idx_t i = 0; i < N; ++i)
            if(i != A_idx)
                isRaised.push_back(A.isRaised(i));

        for(idx_t i = 0; i < M; ++i)
            if(i != B_idx)
                isRaised.push_back(B.isRaised(i));

        POW2_ASSERT(isRaised.size() == ret.raised.size());
        ret.raised = isRaised;

        return ret;
    }

    template<typename T, typename U, typename UU, idx_t N, idx_t M>
    Tensor < typename ENSEM::Promote<typename ENSEM::Promote<T, U>::Type_t, UU>::Type_t, N + M - 2 > contract(const Tensor<T, N> &A,
            const Tensor<U, 2> &metric,
            const Tensor<UU, M> &B,
            const idx_t A_idx,
            const idx_t B_idx)
    {
        if(A.isRaised(A_idx) == B.isRaised(B_idx))
            Tensor<typename ENSEM::Promote<T, U>::Type_t, N> tmp = applyMetric(A, metric, A_idx);

        return contract(A, B, A_idx, B_idx);
    }

    template<typename T, typename U, idx_t N>
    Tensor<typename ENSEM::Promote<T, U>::Type_t , N> operator*(const Tensor<T, N> &lhs, const U rhs)
    {
        Tensor<typename ENSEM::Promote<T, U>::Type_t, N> ret;

        ret.raised = lhs.raised;
        ret.tensor = lhs.tensor * rhs;

        return ret;
    }

    template<typename T, typename U, idx_t N>
    Tensor<typename ENSEM::Promote<T, U>::Type_t , N> operator/(const Tensor<T, N> &lhs, const U rhs)
    {
        Tensor<typename ENSEM::Promote<T, U>::Type_t, N> ret;

        ret.raised = lhs.raised;
        ret.tensor = lhs.tensor / rhs;

        return ret;
    }

    template<typename T, typename U, idx_t N>
    Tensor<typename ENSEM::Promote<T, U>::Type_t , N> operator+(const Tensor<T, N> &lhs, const Tensor<U, N> &rhs)
    {
        Tensor<typename ENSEM::Promote<T, U>::Type_t, N> ret;

        for(idx_t i = 0; i < N; ++i)
            POW2_ASSERT(lhs.isRaised(i) == rhs.isRaised(i));

        ret.raised = lhs.raised;
        ret.tensor = lhs.tensor + rhs.tensor;

        return ret;
    }

    template<typename T, typename U, idx_t N>
    Tensor<typename ENSEM::Promote<T, U>::Type_t , N> operator-(const Tensor<T, N> &lhs, const Tensor<U, N> &rhs)
    {
        return lhs + (-rhs);
    }

    template<typename T, idx_t N>
    std::ostream &operator<<(std::ostream &o, const Tensor<T, N> &t)
    {
        o << "indicies (1 means raised) \n";

        for(idx_t i = 0; i < N; ++i)
            o << "index[" << i << "] = " << t.raised[i] << " ";

        o << "\n";

        o << t.tensor << std::endl;

        return o;
    }

} // close tensor namespace





#endif
