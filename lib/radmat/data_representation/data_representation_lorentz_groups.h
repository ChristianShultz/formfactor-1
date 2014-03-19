#ifndef DATA_REPRESENTATION_LORENTZ_GROUPS_H
#define DATA_REPRESENTATION_LORENTZ_GROUPS_H 

#include "data_representation_primitive_rep.h"


namespace radmat
{
  // base class 
  struct SpherRep_p
  {
    virtual ~SpherRep_p() {}
    virtual int rep_parity() const = 0; 
    virtual int rep_spin() const = 0;
    virtual int rep_eta_p() const = 0; 
  }; 

  // template instantion of parameters 
  template<typename T, int spin, int parity> 
    struct SpherRep
    : public SpherRep_p,
      public LorentzRep_p
    {
      virtual ~SpherRep() {}
      virtual int rep_id() const { return Stringify<T>(); }
      virtual int rep_parity() const { return parity; }
      virtual int rep_spin() const { return spin; }
      virtual int rep_eta_p() const 
      { return (rep_spin() % 2 == 0) ? rep_parity() : -rep_parity();}
    };

  template<typename T> 
  struct SpherRepEmbedType
  : public LorentzRep_p
  { 
    virtual ~SpherRepEmbedType() {}
    virtual std::string rep_id() const { return Stringify<T>(); }
  };

  struct J0p;  
  struct J0m;
  struct J1p;  
  struct J1m;
  struct J2p;  
  struct J2m;
  struct J3p;  
  struct J3m;
  struct J4p;  
  struct J4m;  

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

  struct J0p : public SpherRepEmbedType<J0p> {};  
  struct J0m : public SpherRepEmbedType<J0m> {};
  struct J1p : public SpherRepEmbedType<J1p> {};  
  struct J1m : public SpherRepEmbedType<J1m> {};
  struct J2p : public SpherRepEmbedType<J2p> {};  
  struct J2m : public SpherRepEmbedType<J2m> {};
  struct J3p : public SpherRepEmbedType<J3p> {};  
  struct J3m : public SpherRepEmbedType<J3m> {};
  struct J4p : public SpherRepEmbedType<J4p> {};  
  struct J4m : public SpherRepEmbedType<J4m> {};  

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


  namespace SpherRepresentationFactoryEnv
  {
    bool registerAll(); 
    std::vector<std::string> spher_keys();
    rHandle<SpherRep> callFactory(const std::string &id); 
  }


} 



#endif /* DATA_REPRESENTATION_LORENTZ_GROUPS_H */
