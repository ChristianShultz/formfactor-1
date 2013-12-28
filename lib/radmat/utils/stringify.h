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
      virtual ~StringifyBase() {}
      virtual std::string name() const = 0; 
    }; 

  template<class T>
    struct StringifyType : public StringifyBase
  { 
    ~StringifyType() {}
  };


#define REGISTER_STRINGIFY_TYPE(X)                \
  template<>                                      \
  struct StringifyType<X> : public StringifyBase  \
  {                                               \
    std::string name() const { return #X ;}       \
  };                                              \

  // only specializations may be instatiated
  template<typename T>
    std::string Stringify(void)
    {
      StringifyType<T> f;
      return f.name(); 
    }

}
#endif
