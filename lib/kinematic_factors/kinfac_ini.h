#ifndef KINFAC_INI_H_H_GUARD
#define KINFAC_INI_H_H_GUARD

namespace kinfac
{
  
  /**
     @file kinfac_ini.h 
     @brief a piece of the ini params controlling kinematic factors
   */


  /**
     @brief a piece of the ini params controlling kinematic factors
   */
  struct kinIni
  {
    double massCut;     //! the cut to determine if two states are the same mass 
  };


}
#endif
