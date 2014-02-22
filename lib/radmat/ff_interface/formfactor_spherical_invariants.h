#ifndef FORMFACTOR_CONTINUUM_INVARIANTS_H
#define FORMFACTOR_CONTINUUM_INVARIANTS_H 

#include "formfactor_invariants.h"
#include "radmat/utils/stringify.h"
#include "radmat/utils/handle.h"
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include <string>
#include <sstream> 

namespace radmat
{
  // base class 
  struct SpherRep_p
    : public FFRep_p
  {
    virtual ~SpherRep_p() {}
    virtual std::string rep_id(void) const = 0; 
    virtual int rep_row(void) const = 0; 
    virtual int rep_parity(void) const = 0; 
    virtual int rep_spin(void) const = 0;
    virtual int rep_eta_p(void) const = 0; 
    virtual std::string reg_id(void) const = 0;
  }; 

// template instantion of parameters 
  template<typename T, int spin,  int row, int parity> 
    struct SpherRep
    : public SpherRep_p
    {
      virtual ~SpherRep() {}
      virtual std::string rep_id(void) const { return Stringify<T>(); } 
      virtual int rep_row(void) const { return row; }
      virtual int rep_parity(void) const { return parity; }
      virtual int rep_spin(void) const { return spin; }
      virtual int rep_eta_p(void) const 
      { return (rep_spin() % 2 == 0) ? rep_parity() : -rep_parity();}
      virtual std::string reg_id(void) const
      {
        std::stringstream ss; ss << rep_id() << "_r" << rep_row();
        return ss.str(); 
      }
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

  REGISTER_STRINGIFY_TYPE( J0p ); 
  REGISTER_STRINGIFY_TYPE( J0m ); 
  REGISTER_STRINGIFY_TYPE( J1p ); 
  REGISTER_STRINGIFY_TYPE( J1m ); 
  REGISTER_STRINGIFY_TYPE( J2p ); 
  REGISTER_STRINGIFY_TYPE( J2m ); 
  REGISTER_STRINGIFY_TYPE( J3p ); 
  REGISTER_STRINGIFY_TYPE( J3m ); 

  // one based rows like in adat
  // row = J - h + 1 
  // h = J - row + 1

  // J0 reps
  // <id , spin , row , parity> 
  struct J0p_1 : public SpherRep<J0p,0,1,1> { virtual ~J0p_1() {} };
  struct J0m_1 : public SpherRep<J0m,0,1,-1> { virtual ~J0m_1() {} };

  // J1 reps
  // <id , spin , row , parity> 
  struct J1p_1 : public SpherRep<J1p,1,1,1> { virtual ~J1p_1() {} };
  struct J1m_1 : public SpherRep<J1m,1,1,-1> { virtual ~J1m_1() {} };
  struct J1p_2 : public SpherRep<J1p,1,2,1> { virtual ~J1p_2() {} };
  struct J1m_2 : public SpherRep<J1m,1,2,-1> { virtual ~J1m_2() {} };
  struct J1p_3 : public SpherRep<J1p,1,3,1> { virtual ~J1p_3() {} };
  struct J1m_3 : public SpherRep<J1m,1,3,-1> { virtual ~J1m_3() {} };

  // J2 reps
  // <id , spin , row , parity> 
  struct J2p_1 : public SpherRep<J2p,2,1,1> { virtual ~J2p_1() {} };
  struct J2m_1 : public SpherRep<J2p,2,1,-1> { virtual ~J2m_1() {} };
  struct J2p_2 : public SpherRep<J2p,2,2,1> { virtual ~J2p_2() {} };
  struct J2m_2 : public SpherRep<J2p,2,2,-1> { virtual ~J2m_2() {} };
  struct J2p_3 : public SpherRep<J2p,2,3,1> { virtual ~J2p_3() {} };
  struct J2m_3 : public SpherRep<J2p,2,3,-1> { virtual ~J2m_3() {} };
  struct J2p_4 : public SpherRep<J2p,2,4,1> { virtual ~J2p_4() {} };
  struct J2m_4 : public SpherRep<J2p,2,4,-1> { virtual ~J2m_4() {} };
  struct J2p_5 : public SpherRep<J2p,2,5,1> { virtual ~J2p_5() {} };
  struct J2m_5 : public SpherRep<J2p,2,5,-1> { virtual ~J2m_5() {} };

  // J3 reps
  // <id , spin , row , parity> 
  struct J3p_1 : public SpherRep<J3p,3,1,1> { virtual ~J3p_1() {} };
  struct J3m_1 : public SpherRep<J3p,3,1,-1> { virtual ~J3m_1() {} };
  struct J3p_2 : public SpherRep<J3p,3,2,1> { virtual ~J3p_2() {} };
  struct J3m_2 : public SpherRep<J3p,3,2,-1> { virtual ~J3m_2() {} };
  struct J3p_3 : public SpherRep<J3p,3,3,1> { virtual ~J3p_3() {} };
  struct J3m_3 : public SpherRep<J3p,3,3,-1> { virtual ~J3m_3() {} };
  struct J3p_4 : public SpherRep<J3p,3,4,1> { virtual ~J3p_4() {} };
  struct J3m_4 : public SpherRep<J3p,3,4,-1> { virtual ~J3m_4() {} };
  struct J3p_5 : public SpherRep<J3p,3,5,1> { virtual ~J3p_5() {} };
  struct J3m_5 : public SpherRep<J3p,3,5,-1> { virtual ~J3m_5() {} };
  struct J3p_6 : public SpherRep<J3p,3,6,1> { virtual ~J3p_6() {} };
  struct J3m_6 : public SpherRep<J3p,3,6,-1> { virtual ~J3m_6() {} };
  struct J3p_7 : public SpherRep<J3p,3,7,1> { virtual ~J3p_7() {} };
  struct J3m_7 : public SpherRep<J3p,3,7,-1> { virtual ~J3m_7() {} };

  // factory to hold them all
  typedef Util::SingletonHolder<
    Util::ObjectFactory<SpherRep_p,
			std::string,
			void,
			SpherRep_p* (*)(void),
			Util::StringFactoryError> >
  TheSpherInvariantsFactory;


  namespace SpherInvariantsFactoryEnv
  {
    bool registerAll( void );
    rHandle<SpherRep_p> callFactory(const std::string &id);
    std::vector<std::string> all_keys(void); 
  }
  


} // radmat



#endif /* FORMFACTOR_CONTINUUM_INVARIANTS_H */
