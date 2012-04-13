#ifndef POLARISATION_FACTORY_INVENTORY_H_H_GUARD
#define POLARISATION_FACTORY_INVENTORY_H_H_GUARD

#include "polarisation_factory_keys.h"
#include "tensor/tensorbase.h"
#include <complex>
#include <map>

namespace polarisation
{


  /**
     @file polarisation_factory_inventory
     @brief manages the factory inventory
   */

  /**
     @brief the factory inventory
     @details the inventory is a static map, the keys map to the corresponding polarisation
     tensor which is a pointer to the base type and must be casted back to be of real use

     the implementation is probably not thread safe while writing to the map so we 
     should be sure to only write with one thread 

     the implementation is probably thread safe for reads which suggests that we should 
     simply fill the map once at the beginning of the program in an appropriate place 
     where we have all of the keys
   */
  struct pFacInv
  {
    typedef pFacKey key;
    typedef tensor::TensorImplBase data;
    typedef pFacKeyComp comp;

  public:
    //! do nothing constructor since the data is static
    pFacInv(void);
    //! manage the deletion of the pointers at the end of the entire program
    ~pFacInv(void);

    //! register data to the inventory
    static void pFacInvRegister(const key&, data*e);       // register data to the factory

    //! check if a key is registered
    static bool pFacInvIsRegistered(const key&);           // check if a key is registered

    //! get the data for a key, it must be registered
    static data* pFacInvGetData(const key&);               // get the data for a key

    //! testing method, dumps the inventory into a log file
    static void dumpInventory(void);                       // print the contents of the factory (testing)

  private: // data store
    //! perhaps we should have used an adat singleton?
    static std::map<key,data*,comp> inventory;             // a single static repository to rule them all
  };

}


#endif
