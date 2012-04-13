#ifndef POLARISATION_BASE_H_H_GUARD
#define POLARISATION_BASE_H_H_GUARD

#include "tensor/tensorbase.h"
#include <iostream>

/**
   @file polarisation_factory_keys.h
   @brief contains the definition for the polarization factory keys
  
   @details this header file contains the definition for the polarization factory keys, 
  the key comparator class, overloads on greater (less) than for the keys
  in order to be able to use std map, an overloaded stream operator and 
  an overloaded == (!=) operator
 */


namespace polarisation
{
  using tensor::idx_t;

  /**
     @brief the keys used by the polarisation factory class
   */
  struct pFacKey
  {
    //! construct a key
    pFacKey(tensor::Tensor<double,1> &_pmu, idx_t _J, short _helicity);

    //! assignment operator
    pFacKey& operator=(const pFacKey &o);
    
    //! return pmu[i]
    double& operator[](const idx_t i);

    //! return pmu[i]
    const double& operator[](const idx_t i) const;

    //! key comparison class for std::map
    friend struct pFacKeyComp;

  public: // data store -- public for easy access, its just a container
    tensor::Tensor<double, 1> pmu;   //! the 4 momentum
    idx_t J;                         //! total spin
    short helicity;                  //! maybe its the helicity?
  };
   
  /**
     @brief  key comparison for weak ordering within a std::map
   */
  struct pFacKeyComp
  {
    //! comparison operator
    bool operator() (const pFacKey &lhs, const pFacKey &rhs);
  };

  //! overload of logical < 
  bool operator<(const pFacKey &a, const pFacKey &b);

  //! overload of logical > 
  bool operator>(const pFacKey &a, const pFacKey &b);

  //! overload of logical == 
  bool operator==(const pFacKey &a, const pFacKey &b);

  //! overload of logical != 
  bool operator!=(const pFacKey &a, const pFacKey &b);

  //! overload the stream operator, a key is all public so it doesn't need to  be a friend
  std::ostream& operator<<(std::ostream &o,const pFacKey &k);
  
}

#endif
