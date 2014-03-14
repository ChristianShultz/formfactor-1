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
    virtual std::string rep_id() const = 0;
    virtual std::string rep_g() const = 0; 
    virtual std::string rep_irrep() const = 0; 
    virtual int dim() const = 0; 
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
  template<class G, class IR, int DIM>
  struct CubicRepRest_t
    : public CubicRep_p
  {
    virtual ~CubicRepRest_t() {}; 
    virtual std::string rep_g() const { return Stringify<G>(); }
    virtual std::string rep_irrep() const { return Stringify<IR>(); }
    virtual std::string rep_id() const { return rep_irrep(); }
    virtual int dim() const { return DIM; } 
  };

  
  struct A1Rep_t : public CubicRepRest_t<Oh,A1,1> { virtual ~A1Rep_t() {} };
  struct A2Rep_t : public CubicRepRest_t<Oh,A2,1> { virtual ~A2Rep_t() {} };
  struct T1Rep_t : public CubicRepRest_t<Oh,T1,3> { virtual ~T1Rep_t() {} };
  struct T2Rep_t : public CubicRepRest_t<Oh,T2,3> { virtual ~T2Rep_t() {} };
  struct ERep_t  : public CubicRepRest_t<Oh,E, 2> { virtual ~ERep_t()  {} };



  ////////////////////////////////////////////////////////
  //
  //   FLIGHT IRREPS 
  //

  // lg variants 
  struct E2 : public Rep<E2> { virtual ~E2() {} }; 
  struct B1 : public Rep<B1> { virtual ~B1() {} }; 
  struct B2 : public Rep<B2> { virtual ~B2() {} }; 
  
  REGISTER_STRINGIFY_TYPE(E2);
  REGISTER_STRINGIFY_TYPE(B1);
  REGISTER_STRINGIFY_TYPE(B2);

  // helicities for tokens
  struct H0 : public Rep<H0> { enum {LAMBDA = 0}; virtual ~H0() {} };
  struct H1 : public Rep<H1> { enum {LAMBDA = 1}; virtual ~H1() {} };
  struct H2 : public Rep<H2> { enum {LAMBDA = 2}; virtual ~H2() {} };
  struct H3 : public Rep<H3> { enum {LAMBDA = 3}; virtual ~H3() {} };
  struct H4 : public Rep<H4> { enum {LAMBDA = 4}; virtual ~H4() {} };

  REGISTER_STRINGIFY_TYPE(H0);
  REGISTER_STRINGIFY_TYPE(H1);
  REGISTER_STRINGIFY_TYPE(H2);
  REGISTER_STRINGIFY_TYPE(H3);
  REGISTER_STRINGIFY_TYPE(H4);

  struct CubicRepFlight_p
    : public CubicRep_p
  {
    virtual ~CubicRepFlight_p() {}; 
    virtual std::string rep_h() const = 0; 
    virtual int helicity() const = 0;  
  };


  // rows are 1 based
  template<class HRep, class G, class IR, int DIM>
  struct CubicRepFlight_t
    : public CubicRepFlight_p
  {
    virtual ~CubicRepFlight_t() {}; 
    virtual std::string rep_g() const { return Stringify<G>(); }
    virtual std::string rep_irrep() const { return Stringify<IR>(); }
    virtual std::string rep_h() const { return Stringify<HRep>(); }
    virtual int helicity() const { return HRep::LAMBDA; }
    virtual std::string rep_id() const { return rep_h() + rep_g() + rep_irrep(); }
    virtual int dim() const { return DIM; } 
  };


  // D4 
  struct H0D4A1Rep_t : public CubicRepFlight_t<H0,D4,A1,1> { virtual ~H0D4A1Rep_t() {} };
  struct H0D4A2Rep_t : public CubicRepFlight_t<H0,D4,A2,1> { virtual ~H0D4A2Rep_t() {} };
  struct H1D4E2Rep_t : public CubicRepFlight_t<H1,D4,E2,2> { virtual ~H1D4E2Rep_t() {} };
  struct H2D4B1Rep_t : public CubicRepFlight_t<H2,D4,B1,1> { virtual ~H2D4B1Rep_t() {} };
  struct H2D4B2Rep_t : public CubicRepFlight_t<H2,D4,B2,1> { virtual ~H2D4B2Rep_t() {} };
  struct H3D4E2Rep_t : public CubicRepFlight_t<H3,D4,E2,2> { virtual ~H3D4E2Rep_t() {} };
  struct H4D4A1Rep_t : public CubicRepFlight_t<H4,D4,A1,1> { virtual ~H4D4A1Rep_t() {} };
  struct H4D4A2Rep_t : public CubicRepFlight_t<H4,D4,A2,1> { virtual ~H4D4A2Rep_t() {} };


  // D2 
  struct H0D2A1Rep_t : public CubicRepFlight_t<H0,D2,A1,1> { virtual ~H0D2A1Rep_t() {} };
  struct H0D2A2Rep_t : public CubicRepFlight_t<H0,D2,A2,1> { virtual ~H0D2A2Rep_t() {} };
  struct H1D2B1Rep_t : public CubicRepFlight_t<H1,D2,B1,1> { virtual ~H1D2B1Rep_t() {} };
  struct H1D2B2Rep_t : public CubicRepFlight_t<H1,D2,B2,1> { virtual ~H1D2B2Rep_t() {} };
  struct H2D2A1Rep_t : public CubicRepFlight_t<H2,D2,A1,1> { virtual ~H2D2A1Rep_t() {} };
  struct H2D2A2Rep_t : public CubicRepFlight_t<H2,D2,A2,1> { virtual ~H2D2A2Rep_t() {} };
  struct H3D2B1Rep_t : public CubicRepFlight_t<H3,D2,B1,1> { virtual ~H3D2B1Rep_t() {} };
  struct H3D2B2Rep_t : public CubicRepFlight_t<H3,D2,B2,1> { virtual ~H3D2B2Rep_t() {} };
  struct H4D2A1Rep_t : public CubicRepFlight_t<H4,D2,A1,1> { virtual ~H4D2A1Rep_t() {} };
  struct H4D2A2Rep_t : public CubicRepFlight_t<H4,D2,A2,1> { virtual ~H4D2A2Rep_t() {} };


  // D3 
  struct H0D3A1Rep_t : public CubicRepFlight_t<H0,D3,A1,1> { virtual ~H0D3A1Rep_t() {} };
  struct H0D3A2Rep_t : public CubicRepFlight_t<H0,D3,A2,1> { virtual ~H0D3A2Rep_t() {} };
  struct H1D3E2Rep_t : public CubicRepFlight_t<H1,D3,E2,2> { virtual ~H1D3E2Rep_t() {} };
  struct H2D3E2Rep_t : public CubicRepFlight_t<H2,D3,E2,2> { virtual ~H2D3E2Rep_t() {} };
  struct H3D3A1Rep_t : public CubicRepFlight_t<H3,D3,A1,1> { virtual ~H3D3A1Rep_t() {} };
  struct H3D3A2Rep_t : public CubicRepFlight_t<H3,D3,A2,1> { virtual ~H3D3A2Rep_t() {} };
  struct H4D3E2Rep_t : public CubicRepFlight_t<H4,D3,E2,2> { virtual ~H4D3E2Rep_t() {} };


  // not doing C4 .. too lazy
  

  namespace CubicInvariantsFactoryEnv
  {
    bool registerAll(void); 
    rHandle<CubicRep_p> callFactory(const std::string &s);
    std::vector<std::string> all_keys(void); 
  }

}





#endif /* FORMFACTOR_CUBIC_INVARIANTS_H */
