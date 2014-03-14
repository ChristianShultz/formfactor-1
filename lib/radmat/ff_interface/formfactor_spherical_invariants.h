#ifndef FORMFACTOR_SPHERICAL_INVARIANTS_H
#define FORMFACTOR_SPHERICAL_INVARIANTS_H 

#include "formfactor_invariants.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"
#include <string>
#include <sstream> 

namespace radmat
{
  struct SpherRep_p; 
  REGISTER_STRINGIFY_TYPE(SpherRep_p); 

  // base class 
  struct SpherRep_p
    : public FFRep_p
  {
    virtual ~SpherRep_p() {}
    virtual std::string rep_type(void) const {return Stringify<SpherRep_p>();}
    virtual std::string rep_id(void) const = 0; 
    virtual int rep_parity(void) const = 0; 
    virtual int rep_spin(void) const = 0;
    virtual int rep_eta_p(void) const = 0; 
  }; 

  // template instantion of parameters 
  template<typename T, int spin, int parity> 
    struct SpherRep
    : public SpherRep_p
    {
      virtual ~SpherRep() {}
      virtual std::string rep_id(void) const { return Stringify<T>(); } 
      virtual int rep_parity(void) const { return parity; }
      virtual int rep_spin(void) const { return spin; }
      virtual int rep_eta_p(void) const 
      { return (rep_spin() % 2 == 0) ? rep_parity() : -rep_parity();}
    };

  // spin reps -- this is the rep_id 
  //        this should match the cont_spin class names so we can 
  //        simply smash ids together to form keys for the
  //        lorentzff_formfactor_factory 
  
  struct J0p {};  
  struct J0m {};
  struct J1p {};  
  struct J1m {};
  struct J2p {};  
  struct J2m {};
  struct J3p {};  
  struct J3m {};
  struct J4p {};  
  struct J4m {};  

  REGISTER_STRINGIFY_TYPE( J0p ); 
  REGISTER_STRINGIFY_TYPE( J0m ); 
  REGISTER_STRINGIFY_TYPE( J1p ); 
  REGISTER_STRINGIFY_TYPE( J1m ); 
  REGISTER_STRINGIFY_TYPE( J2p ); 
  REGISTER_STRINGIFY_TYPE( J2m ); 
  REGISTER_STRINGIFY_TYPE( J3p ); 
  REGISTER_STRINGIFY_TYPE( J3m ); 
  REGISTER_STRINGIFY_TYPE( J4p ); 
  REGISTER_STRINGIFY_TYPE( J4m ); 

  struct J0pRep_t : public SpherRep<J0p,0,1>  { virtual ~J0pRep_t() {} };
  struct J0mRep_t : public SpherRep<J0m,0,-1> { virtual ~J0mRep_t() {} };

  struct J1pRep_t : public SpherRep<J1p,1,1>  { virtual ~J1pRep_t() {} };
  struct J1mRep_t : public SpherRep<J1m,1,-1> { virtual ~J1mRep_t() {} };

  struct J2pRep_t : public SpherRep<J2p,2,1>  { virtual ~J2pRep_t() {} };
  struct J2mRep_t : public SpherRep<J2m,2,-1> { virtual ~J2mRep_t() {} };

  struct J3pRep_t : public SpherRep<J3p,3,1>  { virtual ~J3pRep_t() {} };
  struct J3mRep_t : public SpherRep<J3m,3,-1> { virtual ~J3mRep_t() {} };

  struct J4pRep_t : public SpherRep<J4p,4,1>  { virtual ~J4pRep_t() {} };
  struct J4mRep_t : public SpherRep<J4m,4,-1> { virtual ~J4mRep_t() {} };



  namespace SpherInvariantsFactoryEnv
  {
    bool registerAll( void );
    rHandle<SpherRep_p> callFactory(const std::string &id);
    std::vector<std::string> all_keys(void); 
  }
  


} // radmat



#endif /* FORMFACTOR_SPHERICAL_INVARIANTS_H */
