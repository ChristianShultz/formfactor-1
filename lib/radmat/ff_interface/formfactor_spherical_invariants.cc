/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : formfactor_continuum_invariants.cc

* Purpose :

* Creation Date : 22-02-2014

* Last Modified : Sat 22 Feb 2014 04:51:12 PM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "formfactor_spherical_invariants.h"
#include "radmat/utils/pow2assert.h"

namespace FacEnv = radmat::SpherInvariantsFactoryEnv;
typedef radmat::TheSpherInvariantsFactory Factory;


namespace radmat
{

  namespace SpherInvariantsFactoryEnv
  {
    

    bool registered = false; 

    template<class T, class U> 
      T* upCast(void)
      { 
        T *t = new U(); 
        POW2_ASSERT(t);
        return t; 
      }

    template<typename Derived>
      bool
      do_reg(void)
      {
        Derived d; 
        bool reg = Factory::Instance().registerObject( 
            d.reg_id() , upCast<SpherRep_p,Derived> ); 

        if( !!! reg )
          std::cout << __PRETTY_FUNCTION__ 
            << ": reg error for " << d.reg_id() << std::endl;

        return reg; 
      }


    bool registerAll(void)
    {
      bool success = true; 

      if( !!! registered )
      {
        success &= do_reg<J0p_1>(); 
        success &= do_reg<J0m_1>();

        success &= do_reg<J1p_1>(); 
        success &= do_reg<J1m_1>();
        success &= do_reg<J1p_2>(); 
        success &= do_reg<J1m_2>();
        success &= do_reg<J1p_3>(); 
        success &= do_reg<J1m_3>();

        success &= do_reg<J2p_1>(); 
        success &= do_reg<J2m_1>();
        success &= do_reg<J2p_2>(); 
        success &= do_reg<J2m_2>();
        success &= do_reg<J2p_3>(); 
        success &= do_reg<J2m_3>();
        success &= do_reg<J2p_4>(); 
        success &= do_reg<J2m_4>();
        success &= do_reg<J2p_5>(); 
        success &= do_reg<J2m_5>();

        success &= do_reg<J3p_1>(); 
        success &= do_reg<J3m_1>();
        success &= do_reg<J3p_2>(); 
        success &= do_reg<J3m_2>();
        success &= do_reg<J3p_3>(); 
        success &= do_reg<J3m_3>();
        success &= do_reg<J3p_4>(); 
        success &= do_reg<J3m_4>();
        success &= do_reg<J3p_5>(); 
        success &= do_reg<J3m_5>();
        success &= do_reg<J3p_6>(); 
        success &= do_reg<J3m_6>();
        success &= do_reg<J3p_7>(); 
        success &= do_reg<J3m_7>();

        registered = true; 
      }

      return success; 
    }


    // make it blow up if anything goes wrong by wrapping another 
    // call around the factory.createObj method
    rHandle<SpherRep_p > callFactory(const std::string &matElemID)
    {
      SpherRep_p *foo;
      foo = NULL;
      try
      {
        foo = Factory::Instance().createObject(matElemID);
      }
      catch(std::exception &e)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << e.what(); 
        throw e; 
      }
      catch(std::string &s)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << s << std::endl;
        throw s;
      }
      catch(...)
      {
        std::cout << "elem - " << matElemID << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": some error" << std::endl;
        POW2_ASSERT(false); 
      }

      // not a null pointer
      POW2_ASSERT(foo);
      return rHandle<SpherRep_p >(foo);
    }

    // dump all keys in the factory
    std::vector<std::string> 
      all_keys(void)
    {
      return Factory::Instance().keys(); 
    }



  } // SpherInvariantsFactoryEnv

} // radmat

