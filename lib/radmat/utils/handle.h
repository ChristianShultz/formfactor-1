#ifndef HANDLE_INTERFACE_H
#define HANDLE_INTERFACE_H 

#include <stdlib.h>
#include <iostream>

namespace radmat
{

  /*!
   * This is the ADAT handle in a form that doesnt suck 
   */
  template <class T>
    class rHandle
    {
      public:
        //! Initialize pointer with existing pointer
        /*! Requires that the pointer p is a return value of new */
        rHandle(T* p=0) : ptr(p), count(new int(1)) {}

        //! Copy pointer (one more owner)
        rHandle(const rHandle& p) : ptr(p.ptr), count(p.count) 
      {++*count;}

        // dynamic cast via construction 
        template<typename Q>
          rHandle(const rHandle<Q> &p)
          : count(new int(1)) 
          {

            ptr = dynamic_cast<T*>(p.ptr);
            if( ptr == 0x0 ) { 
              std::cerr << "Dynamic cast failed in rHandle::cast()" << std::endl;
              std::cerr << "You are trying to cast to a class you cannot cast to" << std::endl;
              exit(1);
            }
            delete count;

            count = p.count;
            ++*count;
          }


        //! Destructor (delete value if this was the last owner)
        ~rHandle() {dispose();}

        //! Assignment (unshare old and share new value)
        rHandle& operator=(const rHandle& p) 
        {
          if (this != &p) 
          {
            dispose();
            ptr = p.ptr;
            count = p.count;
            ++*count;
          }
          return *this;
        }


        //! The cast function requires all rHandles<Q> to be friends of rHandle<T>
        template<typename Q> friend class rHandle;

        //! Access the value to which the pointer refers
        T& operator*() const {return *ptr;}
        T* operator->() const {return ptr;}
        T* get_ptr() const {return ptr;}

      private:
        void dispose() 
        {
          if (--*count == 0) 
          {
            delete count;
            delete ptr;
          }
        }

      private:
        T* ptr;        // pointer to the value
        mutable int* count;    // shared number of owners
    };

  //  // this is a nasty hack to get temporary access to the 
  //  // raw pointer that lives inside an adat handle, 
  //  // as always one should be very careful if using it
  //  template<typename Derived, typename Base> 
  //    const rHandle<Derived>
  //    dynamic_cast_handle(const rHandle<Base> &h)
  //    {
  //      return h.cast<Derived>();
  //    }


} // radmat 


#endif /* HANDLE_INTERFACE_H */
