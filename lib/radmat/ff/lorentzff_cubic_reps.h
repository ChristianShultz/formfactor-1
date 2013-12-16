#ifndef LORENTZFF_CUBIC_REPS_H
#define LORENTZFF_CUBIC_REPS_H 

#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"
#include "hadron/irrep_util.h"
#include "io/adat_xmlio.h"

namespace radmat
{

  struct Rep_p
  {
    virtual ~Rep_p() {}
    virtual std::string id(void) const = 0; 
  };  


  template<typename T> 
    struct Rep
    : public Rep_p
    {
      virtual ~Rep() {}
      virtual std::string id(void) const {return Stringify<T>();}
    };

  // objects carry type info 
  struct Oh : public Rep<Oh> { virtual ~Oh() {} };
  struct D2 : public Rep<D2> { virtual ~D2() {} };
  struct D3 : public Rep<D3> { virtual ~D3() {} };
  struct D4 : public Rep<D4> { virtual ~D4() {} };

  REGISTER_STRINGIFY_TYPE(Oh);
  REGISTER_STRINGIFY_TYPE(D2);
  REGISTER_STRINGIFY_TYPE(D3);
  REGISTER_STRINGIFY_TYPE(D4);

  struct RepHandler
  {

    rHandle<Rep_p> gen_rep(const ADATXML::Array<int> &p) const
    {
      std::string LG = Hadron::generateLittleGroup(p); 
      if(LG == "Oh")
        return rHandle<Rep_p>(new Oh() ); 
      else if (LG == "D2")
        return rHandle<Rep_p>(new D2() ); 
      else if (LG == "D3")
        return rHandle<Rep_p>(new D3() ); 
      else if (LG == "D4")
        return rHandle<Rep_p>(new D4() ); 
      else
      {
        std::cout << "LG " << LG << " not supported" << std::endl;
        exit(1); 
      }
      exit(1); 
    }

  };

  struct RepPair
  {
    RepPair(const ADATXML::Array<int> &left, const ADATXML::Array<int> &right)
      : l(left) , r(right)
    { 
      RepHandler R; 
      lefty = R.gen_rep(l);
      righty = R.gen_rep(r); 
    }

    ADATXML::Array<int> l;
    ADATXML::Array<int> r; 
    rHandle<Rep_p> lefty;
    rHandle<Rep_p> righty; 
  }; 

}
#endif /* LORENTZFF_CUBIC_REPS_H */
