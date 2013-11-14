#ifndef HANDLE_INTERFACE_H
#define HANDLE_INTERFACE_H 

#include "adat/handle.h"


namespace radmat
{
  
  // this is a nasty hack to get temporary access to the 
  // raw pointer that lives inside an adat handle, 
  // as always one should be very careful if using it
  template<typename Derived, typename Base> 
    const Derived *
    dynamic_cast_handle(const ADAT::Handle<Base> &h)
    {
      const Derived* t = dynamic_cast<const Derived*>(&*h);
    }


} // radmat 


#endif /* HANDLE_INTERFACE_H */
