#ifndef REDSTAR_UNIMPROVED_VECTOR_CURRENT_H
#define REDSTAR_UNIMPROVED_VECTOR_CURRENT_H 

#include "redstar_abstract_block_base.h"
#include "radmat/utils/stringify.h"
#include "io/adat_xmlio.h"
#include <utility>

namespace radmat
{

  // a photon fragment 
  struct RedstarUnimprovedVectorCurrentPFrag
  {
    double coeff; 
    std::string name; 
  };

  //
  //
  // A BASIC gamma^mu VECTOR CURRENT INPUT
  //
  //

  struct RedstarUnimprovedVectorCurrentInput;
  REGISTER_STRINGIFY_TYPE(RedstarUnimprovedVectorCurrentInput);

  struct RedstarUnimprovedVectorCurrentInput
    : public AbsRedstarInput_t 
  {
    typedef RedstarUnimprovedVectorCurrentPFrag pfrag;
    typedef std::vector<pfrag>::const_iterator const_iterator; 

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarUnimprovedVectorCurrentInput>()); 
    }

    virtual std::string write(void) const; 

    virtual rHandle<AbsRedstarInput_t>  clone(void) const 
    {
      return rHandle<AbsRedstarInput_t>(
          new RedstarUnimprovedVectorCurrentInput(*this) ); 
    }

    const_iterator begin(void) const {return photons.begin();}
    const_iterator end(void) const {return photons.end();}

    int lorentz;                    // 1,2,3,4 -- euclidean so time is 4 
    ADATXML::Array<int> mom;        // momentum
    std::vector<pfrag> photons;     // some photon pieces  
    bool creation_op;               // is it a creation op
    bool smearedP;                  // is it smeared
    int t_slice;                    // where does it live
  };



  //
  //
  // A BASIC gamma^mu VECTOR CURRENT FUNCTOR
  //
  //

  struct RedstarUnimprovedVectorCurrentBlock;
  REGISTER_STRINGIFY_TYPE(RedstarUnimprovedVectorCurrentBlock);

  struct RedstarUnimprovedVectorCurrentBlock
    : public AbsRedstarBlock_t
  {
    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarUnimprovedVectorCurrentBlock>()); 
    }

    virtual EnsemRedstarBlock
      operator()(const AbsRedstarInput_t *ptr2derived) const ; 

  };


  //
  //
  // A BASIC gamma^mu VECTOR CURRENT XML INTERFACE
  //
  //


  struct RedstarUnimprovedVectorCurrentXML;
  REGISTER_STRINGIFY_TYPE(RedstarUnimprovedVectorCurrentXML);


  struct RedstarUnimprovedVectorCurrentXML 
    : public AbsRedstarXMLInterface_t 
  {
    typedef RedstarUnimprovedVectorCurrentPFrag pfrag;

    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path); 

    virtual std::string type(void) const
    {
      return std::string(Stringify<RedstarUnimprovedVectorCurrentXML>());
    }

    virtual std::string write(void) const; 
    virtual void write(ADATXML::XMLWriter &xml, const std::string &path) const; 

    virtual int timeslice(void) const { return t_slice; }

    struct insertion
    {
      bool active; 
      bool creation_op; 
      bool smearedP; 
      ADATXML::Array<pfrag> photons;  
    };

    int pmin; 
    int pmax; 
    int t_slice; 
    insertion time; 
    insertion space; 
  };

  // some xml trash
  void read(ADATXML::XMLReader &xml, 
      const std::string &path,
      RedstarUnimprovedVectorCurrentPFrag &p);

  void write(ADATXML::XMLWriter &xml, 
      const std::string &path,
      const RedstarUnimprovedVectorCurrentPFrag &p);

  void read(ADATXML::XMLReader &xml, 
      const std::string &path,
      RedstarUnimprovedVectorCurrentXML::insertion &i);

  void write(ADATXML::XMLWriter &xml,
      const std::string &path, 
      const RedstarUnimprovedVectorCurrentXML::insertion &i);
} // radmat


#endif /* REDSTAR_UNIMPROVED_VECTOR_CURRENT_H */
