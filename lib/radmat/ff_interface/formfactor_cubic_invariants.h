#ifndef FORMFACTOR_CUBIC_INVARIANTS_H
#define FORMFACTOR_CUBIC_INVARIANTS_H 


#include "radmat/ff/lorentzff_cubic_reps.h"
#include "formfactor_invariants.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/stringify.h"
#include <sstream>
#include <string> 

// make classes reflecting the properties of the irreps 

namespace radmat
{
  
  struct CubicRep_p;
  REGISTER_STRINGIFY_TYPE(CubicRep_p); 

  struct CubicRep_p
    : public FFRep_p
  {
    virtual ~CubicRep_p() {}
    virtual std::string rep_type() const { return Stringify<CubicRep_p>(); }
    virtual std::string rep_id() const { return rep_token(); }
    virtual std::string rep_g() const = 0; 
    virtual std::string rep_irrep() const = 0; 
    virtual std::string rep_token() const = 0; 
    virtual int dim() const = 0; 
    virtual int rep_row() const = 0; 

    virtual std::string reg_id() const = 0; 
    std::string itos(const int i) const { std::stringstream ss; ss << i; return ss.str(); }
  };

  
  ////////////////////////////////////////////////////////
  //
  //   REST IRREPS 
  //


  // all the irreps -- Oh
  struct A1 : public Rep<A1> { virtual ~A1() {} };
  struct A2 : public Rep<A2> { virtual ~A2() {} };
  struct T1 : public Rep<T1> { virtual ~T1() {} };
  struct T2 : public Rep<T2> { virtual ~T2() {} };
  struct E : public Rep<E> { virtual ~E() {} };

  REGISTER_STRINGIFY_TYPE(A1);
  REGISTER_STRINGIFY_TYPE(A2);
  REGISTER_STRINGIFY_TYPE(T1);
  REGISTER_STRINGIFY_TYPE(T2);
  REGISTER_STRINGIFY_TYPE(E); 

  // rows are 1 based
  template<class G, class IR, int DIM, int ROW>
  struct CubicRepRest_t
    : public CubicRep_p
  {
    virtual ~CubicRepRest_t() {}; 
    virtual std::string rep_g() const { return Stringify<G>(); }
    virtual std::string rep_irrep() const { return Stringify<IR>(); }
    virtual std::string rep_token() const { return rep_irrep(); }
    virtual int dim() const { return DIM; } 
    virtual int rep_row() const { return ROW; }
    virtual std::string reg_id() const 
    {
      return rep_token() + "_r" + CubicRep_p::itos( rep_row() ) ; 
    }
  };

  
  struct A1Rep_r1 : public CubicRepRest_t<Oh,A1,1,1> { virtual ~A1Rep_r1() {} };
  struct A2Rep_r1 : public CubicRepRest_t<Oh,A2,1,1> { virtual ~A2Rep_r1() {} };
  struct T1Rep_r1 : public CubicRepRest_t<Oh,T1,1,1> { virtual ~T1Rep_r1() {} };
  struct T1Rep_r2 : public CubicRepRest_t<Oh,T1,1,2> { virtual ~T1Rep_r2() {} };
  struct T1Rep_r3 : public CubicRepRest_t<Oh,T1,1,3> { virtual ~T1Rep_r3() {} };
  struct T2Rep_r1 : public CubicRepRest_t<Oh,T2,1,1> { virtual ~T2Rep_r1() {} };
  struct T2Rep_r2 : public CubicRepRest_t<Oh,T2,1,2> { virtual ~T2Rep_r2() {} };
  struct T2Rep_r3 : public CubicRepRest_t<Oh,T2,1,3> { virtual ~T2Rep_r3() {} };
  struct ERep_r1 : public CubicRepRest_t<Oh,E,1,1> { virtual ~ERep_r1() {} };
  struct ERep_r2 : public CubicRepRest_t<Oh,E,1,2> { virtual ~ERep_r2() {} };



  ////////////////////////////////////////////////////////
  //
  //   FLIGHT IRREPS 
  //

  // lg variants 
  struct E1 : public Rep<E1> { virtual ~E1() {} }; 
  struct E2 : public Rep<E2> { virtual ~E2() {} }; 
  struct E3 : public Rep<E3> { virtual ~E3() {} }; 
  struct B1 : public Rep<B1> { virtual ~B1() {} }; 
  struct B2 : public Rep<B2> { virtual ~B2() {} }; 
  
  REGISTER_STRINGIFY_TYPE(E1);
  REGISTER_STRINGIFY_TYPE(E2);
  REGISTER_STRINGIFY_TYPE(E3);
  REGISTER_STRINGIFY_TYPE(B1);
  REGISTER_STRINGIFY_TYPE(B2);

  // helicities for tokens
  struct H0 : public Rep<H0> { virtual ~H0() {} };
  struct H1 : public Rep<H1> { virtual ~H1() {} };
  struct H2 : public Rep<H2> { virtual ~H2() {} };
  struct H3 : public Rep<H3> { virtual ~H3() {} };
  struct H4 : public Rep<H4> { virtual ~H4() {} };

  REGISTER_STRINGIFY_TYPE(H0);
  REGISTER_STRINGIFY_TYPE(H1);
  REGISTER_STRINGIFY_TYPE(H2);
  REGISTER_STRINGIFY_TYPE(H3);
  REGISTER_STRINGIFY_TYPE(H4);


  // rows are 1 based
  template<class HRep, class G, class IR, int DIM, int ROW>
  struct CubicRepFlight_t
    : public CubicRep_p
  {
    virtual ~CubicRepFlight_t() {}; 
    virtual std::string rep_g() const { return Stringify<G>(); }
    virtual std::string rep_irrep() const { return Stringify<IR>(); }
    virtual std::string rep_h() const { return Stringify<HRep>(); }
    virtual std::string rep_token() const { return rep_h() + rep_g() + rep_irrep(); }
    virtual int dim() const { return DIM; } 
    virtual int rep_row() const { return ROW; }
    virtual std::string reg_id() const 
    {
      return rep_token() + "_r" + CubicRep_p::itos( rep_row() ) ; 
    }
  };


  // D4 
  struct H0D4A1Rep_r1 : public CubicRepFlight_t<H0,D4,A1,1,1> { virtual ~H0D4A1Rep_r1() {} };
  struct H0D4A2Rep_r1 : public CubicRepFlight_t<H0,D4,A2,1,1> { virtual ~H0D4A2Rep_r1() {} };

  struct H1D4E2Rep_r1 : public CubicRepFlight_t<H1,D4,E2,2,1>  { virtual ~H1D4E2Rep_r1() {} };
  struct H1D4E2Rep_r2 : public CubicRepFlight_t<H1,D4,E2,2,2>  { virtual ~H1D4E2Rep_r2() {} };

  struct H2D4B1Rep_r1 : public CubicRepFlight_t<H2,D4,B1,1,1> { virtual ~H2D4B1Rep_r1() {} };
  struct H2D4B2Rep_r1 : public CubicRepFlight_t<H2,D4,B2,1,1> { virtual ~H2D4B2Rep_r1() {} };

  struct H3D4E2Rep_r1 : public CubicRepFlight_t<H3,D4,E2,2,1> { virtual ~H3D4E2Rep_r1() {} };
  struct H3D4E2Rep_r2 : public CubicRepFlight_t<H3,D4,E2,2,2> { virtual ~H3D4E2Rep_r2() {} };

  struct H4D4A1Rep_r1 : public CubicRepFlight_t<H4,D4,A1,1,1> { virtual ~H4D4A1Rep_r1() {} };
  struct H4D4A2Rep_r1 : public CubicRepFlight_t<H4,D4,A2,1,1> { virtual ~H4D4A2Rep_r1() {} };


  // D2 
  struct H0D2A1Rep_r1 : public CubicRepFlight_t<H0,D2,A1,1,1> { virtual ~H0D2A1Rep_r1() {} };
  struct H0D2A2Rep_r1 : public CubicRepFlight_t<H0,D2,A2,1,1> { virtual ~H0D2A2Rep_r1() {} };

  struct H1D2B1Rep_r1 : public CubicRepFlight_t<H1,D2,B1,1,1> { virtual ~H1D2B1Rep_r1() {} };
  struct H1D2B2Rep_r1 : public CubicRepFlight_t<H1,D2,B2,1,1> { virtual ~H1D2B2Rep_r1() {} };

  struct H2D2A1Rep_r1 : public CubicRepFlight_t<H2,D2,A1,1,1> { virtual ~H2D2A1Rep_r1() {} };
  struct H2D2A2Rep_r1 : public CubicRepFlight_t<H2,D2,A2,1,1> { virtual ~H2D2A2Rep_r1() {} };

  struct H3D2B1Rep_r1 : public CubicRepFlight_t<H3,D2,B1,1,1> { virtual ~H3D2B1Rep_r1() {} };
  struct H3D2B2Rep_r1 : public CubicRepFlight_t<H3,D2,B2,1,1> { virtual ~H3D2B2Rep_r1() {} };

  struct H4D2A1Rep_r1 : public CubicRepFlight_t<H4,D2,A1,1,1> { virtual ~H4D2A1Rep_r1() {} };
  struct H4D2A2Rep_r1 : public CubicRepFlight_t<H4,D2,A2,1,1> { virtual ~H4D2A2Rep_r1() {} };


  // D3 
  struct H0D3A1Rep_r1 : public CubicRepFlight_t<H0,D3,A1,1,1> { virtual ~H0D3A1Rep_r1() {} };
  struct H0D3A2Rep_r1 : public CubicRepFlight_t<H0,D3,A2,1,1> { virtual ~H0D3A2Rep_r1() {} };

  struct H1D3E2Rep_r1 : public CubicRepFlight_t<H1,D3,E2,2,1> { virtual ~H1D3E2Rep_r1() {} };
  struct H1D3E2Rep_r2 : public CubicRepFlight_t<H1,D3,E2,2,2> { virtual ~H1D3E2Rep_r2() {} };

  struct H2D3E2Rep_r1 : public CubicRepFlight_t<H2,D3,E2,2,1> { virtual ~H2D3E2Rep_r1() {} };
  struct H2D3E2Rep_r2 : public CubicRepFlight_t<H2,D3,E2,2,2> { virtual ~H2D3E2Rep_r2() {} };

  struct H3D3A1Rep_r1 : public CubicRepFlight_t<H3,D3,A1,1,1> { virtual ~H3D3A1Rep_r1() {} };
  struct H3D3A2Rep_r1 : public CubicRepFlight_t<H3,D3,A2,1,1> { virtual ~H3D3A2Rep_r1() {} };

  struct H4D3E2Rep_r1 : public CubicRepFlight_t<H4,D3,E2,2,1> { virtual ~H4D3E2Rep_r1() {} };
  struct H4D3E2Rep_r2 : public CubicRepFlight_t<H4,D3,E2,2,2> { virtual ~H4D3E2Rep_r2() {} };


  // not doing C4 .. too lazy
  

  namespace CubicInvariantsFactoryEnv
  {
    bool registerAll(void); 
    rHandle<CubicRep_p> callFactory(const std::string &s);
    std::vector<std::string> all_keys(void); 
  }

}





#endif /* FORMFACTOR_CUBIC_INVARIANTS_H */
