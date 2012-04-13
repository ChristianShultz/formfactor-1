// polarisation_factory_inventory.cc -
//
// Thursday, March 29 2012
//

#include "polarisation_factory_inventory.h"
#include "polarisation_factory_keys.h"
#include <fstream>
#include <map>


using namespace polarisation;

typedef pFacKey key;
typedef tensor::TensorImplBase data;
typedef pFacKeyComp comp;

std::map<key,data*,comp> pFacInv::inventory;

pFacInv::~pFacInv(void)
{
  std::map<key,data*,comp>::iterator it;

    for (it = inventory.begin(); it != inventory.end(); ++it)  
      delete it->second;    

  inventory.clear();
}

void pFacInv::pFacInvRegister(const key &k,data *d)
{
  // this should never happen
  if(inventory.find(k) != inventory.end())
    std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << "\n"
	      << "Warning, re-registering data" << std::endl;

  POW2_ASSERT(d);

  inventory[k] = d;
}

bool pFacInv::pFacInvIsRegistered(const key &k)
{
  return (inventory.find(k) != inventory.end());
}

data* pFacInv::pFacInvGetData(const key &k)
{
  POW2_ASSERT(pFacInvIsRegistered(k));
  
  data * foo = inventory.find(k)->second;

  POW2_ASSERT(foo);
  
  return foo;
}

void pFacInv::dumpInventory(void) 
{
  std::ofstream out;
  out.open("pFacInvDump.log", std::ios::out | std::ios::app);

  std::map<key,data*,comp>::const_iterator it;

  unsigned int ct(0);

  for (it = inventory.begin(); it != inventory.end(); ++it)  
    {
      ++ct;
      out << "data = " << it->second << "\n";
      out << "key = " << it->first;
    }

  out << "the factory contained " << ct << "tensors" << std::endl;
}
