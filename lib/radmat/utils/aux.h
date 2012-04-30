#ifndef AUX_H_H_GUARD
#define AUX_H_H_GUARD

#include "pow2assert.h"

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


}

#endif
