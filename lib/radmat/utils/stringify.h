#ifndef STRINGIFY_H_H_GUARD
#define STRINGIFY_H_H_GUARD

#include <complex>
#include <string>

/**
  @file stringify.h

  @brief this template is used for testing and prints the type that was tested
  @details ie cout << Stringify<double>() << endl;
  will print double to the terminal, its stupid but its nice for the testing
  */

namespace radmat
{


    struct StringifyBase
    {
      StringifyBase() {}
      virtual std::string name() const = 0; 
    }; 

  template<class T>
    struct StringifyTemp : public StringifyBase
  { };


#define REGISTER_STRINGIFY_TYPE(X)                \
  template<>                                      \
  struct StringifyTemp<X> : public StringifyBase  \
  {                                               \
    std::string name() const { return #X ;}       \
  };                                              \

  // only specializations may be instatiated
  template<typename T>
    std::string Stringify(void)
    {
      StringifyTemp<T> f;
      return f.name(); 
    }


//  template<typename T>
//    struct TypeName
//    {
//      static const char *name;
//    };
//
//  template<typename T>
//    const char *TypeName<T>::name = "unknown";
//
//  template<typename T>
//   inline const char *Stringify(void)
//    {
//      return TypeName<T>::name;
//    }
//
//  // macro template specialization expansion
//#define REGISTER_STRINGIFY_TYPE(X) template<>	\
//  const char *TypeName<X>::name = #X
//
//  // types must be registered before use otherwise it will default to
//  // the unknown type
//  REGISTER_STRINGIFY_TYPE(int);
//  REGISTER_STRINGIFY_TYPE(float);
//  REGISTER_STRINGIFY_TYPE(double);
//  REGISTER_STRINGIFY_TYPE(double);
//  REGISTER_STRINGIFY_TYPE(std::complex<double>);

}
#endif
