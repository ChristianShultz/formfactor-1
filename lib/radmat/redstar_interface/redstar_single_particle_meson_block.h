#ifndef REDSTAR_SINGLE_PARTICLE_MESON_BLOCK_H
#define REDSTAR_SINGLE_PARTICLE_MESON_BLOCK_H 


#include "redstar_abstract_block_base.h"
#include "radmat/utils/stringify.h"

namespace radmat
{

  // !!this is a continuum object!!
  // specify an operator corresponding to a single type, spin, helicity, and momentum
  //    -- may be a sum over lattice values 



  // 
  //
  // SINGLE PARTICLE MESON INPUT LIST OBJECT
  //
  //

  struct RedstarSingleParticleMesonInput_string {};
  REGISTER_STRINGIFY_TYPE(RedstarSingleParticleMesonInput_string); 

  struct RedstarSingleParticleMesonInput
    : public AbsRedstarInput_t
  {
    virtual std::string type(void) const 
    {
      return std::string(Stringify<RedstarSingleParticleMesonInput_string>()); 
    }

    virtual std::string write(void) const; 

    int J;                            // continuum spin
    int H;                            // helicity 
    bool parity;                      // true is positive
    ADATXML::Array<int> mom;          // momentum 
    int twoI_z;                       // twice Isospin 
    std::string name;                 // particle stem 
    bool creation_op;                 // is it a creation operaor
    bool smearedP;                    // is it a smeared operator
    bool isProjected;                 // is it a projected operator
    int t_slice;                      // t_slice
  };


  // 
  //
  // SINGLE PARTICLE MESON FUNCTOR
  //
  //

  struct RedstarSingleParticleMesonBlock_string {};
  REGISTER_STRINGIFY_TYPE(RedstarSingleParticleMesonBlock_string); 

  struct RedstarSingleParticleMesonBlock
    : public AbsRedstarBlock_t
  {
    virtual std::string type(void) const 
    {
      return std::string(Stringify<RedstarSingleParticleMesonBlock_string>()); 
    }

    // arg is polymorphic dynamic cast to derived type RedstarSingleParticleMesonInput
    virtual EnsemRedstarBlock 
      operator()(const AbsRedstarInput_t * ptr2derived) const ; 
  };





  // 
  //
  // SINGLE PARTICLE MESON XML INTERFACE
  //
  //


  struct RedstarSingleParticleMesonXML_string {}; 
  REGISTER_STRINGIFY_TYPE(RedstarSingleParticleMesonXML_string); 


  struct RedstarSingleParticleMesonXML
    : public AbsRedstarXMLInterface_t
  {

    virtual std::string type(void) const 
    {
      return std::string(Stringify<RedstarSingleParticleMesonXML_string>()); 
    }

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string write(void) const;  
    virtual void write(ADATXML::XMLWriter &xml, 
        const std::string &path); 

    int J;                                              // continuum spin
    ADATXML::Array<int> H;                              // helicity 
    bool parity;                                        // true is positive
    bool fill_star;                                     // consider all rotations?
    ADATXML::Array< ADATXML::Array<int> > mom;          // momentum 
    int twoI_z;                                         // twice Isospin 
    std::string name;                                   // particle stem 
    bool creation_op;                                   // is it a creation operaor
    bool smearedP;                                      // is it a smeared operator
    bool isProjected;                                   // is it a projected operator
    int t_slice;                                        // t_slice
  }; 










} // radmat


#endif /* REDSTAR_SINGLE_PARTICLE_MESON_BLOCK_H */
