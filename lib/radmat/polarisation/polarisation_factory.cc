// polarisation_factory.cc -
//
// Saturday, March 31 2012
//

#include "polarisation_factory.h"
#include "radmat/tensor/tensorbase.h"
#include "radmat/utils/pow2assert.h"
#include "polarisation_factory_inventory.h"
#include <iostream>

namespace radmat 
{
  typedef pFacKey key;
  typedef TensorImplBase data;

  pFac::pFac(void)
  {  }


  data* pFac::get(const key &k)
  {
    if(pFacInv::pFacInvIsRegistered(k))
      return pFacInv::pFacInvGetData(k);

    data* foo;

    // a nasty switch since we can't use a variable as a template parameter.
    // fortunately the variable parameter is a prtimitive type and lines up nicely
    // with the case index
    switch(k.J)
      {
      case 0:
	POW2_ASSERT(false);      // this doesn't make sense
	break;
      case 1:
	{
	  Coupler<1> cc(k,this);
	  foo = cc();
	  pFacInv::pFacInvRegister(k,foo);
	  break;
	}
      case 2:
	{
	  Coupler<2> cc(k,this);
	  foo = cc();
	  pFacInv::pFacInvRegister(k,foo);
	  break;
	}
      case 3:
	{
	  Coupler<3> cc(k,this);
	  foo = cc();
	  pFacInv::pFacInvRegister(k,foo);
	  break;
	}
      case 4:
	{
	  Coupler<4> cc(k,this);
	  foo = cc();
	  pFacInv::pFacInvRegister(k,foo);
	  break;
	}
      case 5:
	{
	  Coupler<5> cc(k,this);
	  foo = cc();
	  pFacInv::pFacInvRegister(k,foo);
	  break;
	}
      default:
	POW2_ASSERT(false);     // I don't think we need to go past 4..
      }                         

    return foo;
  }

  void pFac::dumpFactory(void) const
  {
    pFacInv::dumpInventory();
  }

}
