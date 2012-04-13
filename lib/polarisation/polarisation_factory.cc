// polarisation_factory.cc -
//
// Saturday, March 31 2012
//

#include "polarisation_factory.h"
#include "tensor/tensorbase.h"
#include "utils/pow2assert.h"
#include "polarisation_factory_inventory.h"
#include <iostream>

using namespace polarisation;

typedef pFacKey key;
typedef tensor::TensorImplBase data;

pFac::pFac(void)
{  }


data* pFac::get(const key &k)
{
  if(polarisation::pFacInv::pFacInvIsRegistered(k))
    return polarisation::pFacInv::pFacInvGetData(k);

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
	polarisation::pFacInv::pFacInvRegister(k,foo);
	break;
      }
    case 2:
      {
	Coupler<2> cc(k,this);
	foo = cc();
	polarisation::pFacInv::pFacInvRegister(k,foo);
	break;
      }
    case 3:
      {
	Coupler<3> cc(k,this);
	foo = cc();
	polarisation::pFacInv::pFacInvRegister(k,foo);
	break;
      }
    case 4:
      {
	Coupler<4> cc(k,this);
	foo = cc();
	polarisation::pFacInv::pFacInvRegister(k,foo);
	break;
      }
    case 5:
      {
	Coupler<5> cc(k,this);
	foo = cc();
	polarisation::pFacInv::pFacInvRegister(k,foo);
	break;
      }
    default:
      POW2_ASSERT(false);     // if you need more for some perverse reason then just add up to the corresponding case, 
    }                         // I suggest using perl to write the code for you if J is large. 

  return foo;
}

void pFac::dumpFactory(void) const
{
  polarisation::pFacInv::dumpInventory();
}
