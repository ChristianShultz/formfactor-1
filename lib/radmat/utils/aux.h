#ifndef AUX_H_H_GUARD
#define AUX_H_H_GUARD

#include "pow2assert.h"
#include "type_computations.h"

namespace radmat
{

  template<class T, class U >
  T downcastAndDeref(const U *ptr)
  {
    const T *p = dynamic_cast<const T*>(ptr);  // cast
    POW2_ASSERT(p);                            // check cast
    T ret = *p;                                // make a T
    return ret; 
  }

  // stl complex is a pain.. , use this to get some asymetric complex operatons done

  template<class T,class U>
  inline typename Promote<T, U>::Type_t convert_stl_type(const T &t)
  {
    return typename Promote<T, U>::Type_t(t);
  }

  template<typename T, typename U, class BinaryOperator>
  inline typename Promote<T, U>::Type_t binary_op(const T &t, const U &u, BinaryOperator b_op)
  {
    return b_op(convert_stl_type<T,U>(t),convert_stl_type<U,T>(u));
  }

}

#endif
