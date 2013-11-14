#ifndef REDSTAR_ABSTRACT_BLOCK_BASE_H
#define REDSTAR_ABSTRACT_BLOCK_BASE_H 

#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/obj_expr_t.h"
#include "ensem/ensem.h"
#include <string>

namespace radmat
{

  // the list object expressions that we will be playing with, these represent 
  // actual sums over correlation functions that we will be computing
  typedef ListObjExpr_t<ENSEM::Complex,
          Hadron::KeyHadronNPartNPtCorr_t::NPoint_t> EnsemRedstarBlock;

  // an abstract polymorphic base to cast from 
  //    this is a container that gets eaten by the block types, 
  //    the actual impls have the data we are going to need
  struct 
    AbsRedstarInput_t
    {
      virtual std::string type(void) const = 0; 
      virtual std::string write(void) const = 0; 
      virtual ~AbsRedstarInput_t() {}
      virtual AbsRedstarInput_t * clone(void) const = 0; 
    }; 

  // an abstract polymorphic base to cast from
  //   this is an object functor that knows how to make EnsemRedstarBlocks 
  //   provided some abstract input type that corresponds to some 
  //   type of quark current object that redstar knows about
  
  template<class Derived> 
  struct
    AbsRedstarBlock_t
    {
      virtual ~AbsRedstarBlock_t() {}
      virtual AbsRedstarBlock_t * clone(void) const = 0; 
      virtual std::string type(void) const = 0; 
      virtual EnsemRedstarBlock operator()(const AbsRedstarInput_t * ) const = 0; 
      virtual EnsemRedstarBlock operator()(const ADAT::Handle<AbsRedstarInput_t> &o) const
      {
        const AbsRedstarInput_t * ptr = &*o; 
        return static_cast<Derived*>(this)->operator()(ptr); 
      }
    };


  template<class DerivedBlock_t>
  struct AbsRedstarXMLInterface_t
  {
    typedef std::vector<const AbsRedstarInput_t*>::const_iterator const_iterator; 

    AbsRedstarXMLInterface_t() {}
    AbsRedstarXMLInterface_t& operator=(const AbsRedstarXMLInterface_t &); // not impl

    virtual AbsRedstarXMLInterface_t * clone(void) const = 0; 

    // teaches us how to read a given xml type and populate 
    // the input list / allocate memory for the functor & list
    virtual void read(ADATXML::XMLReader &xml, 
        const std::string &path) = 0;     
    virtual std::string write(void) const = 0; 
    virtual std::string type(void) const = 0; 
    virtual void write(ADATXML::XMLWriter &xml, const std::string &path) const = 0; 

    // all npoints have some common info -- this happens to be the most useful 
    virtual int timeslice(void) const = 0;   

    // some stl like overloads
    virtual const_iterator begin(void) const {return inputList.begin();}
    virtual const_iterator end(void) const {return inputList.end();}
    virtual bool empty(void) const {return inputList.empty();}

    // a nice way of storing information and then getting out the continuum version
    virtual EnsemRedstarBlock operator()(const AbsRedstarInput_t *inp) { (*objFunctorPtr)(inp); }

    typename ADAT::Handle< AbsRedstarBlock_t<DerivedBlock_t> > objFunctorPtr; 
    std::vector< ADAT::Handle< AbsRedstarInput_t > > inputList; 
  }; 



} // radmat

#endif /* REDSTAR_ABSTRACT_BLOCK_BASE_H */
