#ifndef ABSTRACT_FUNCTION_PARAMETER_NAMES_H
#define ABSTRACT_FUNCTION_PARAMETER_NAMES_H 


#include <string>
#include <vector>
#include <sstream>
#include "abstract_function_exit.h"
#include "abstract_function_parameter_container.h"

namespace jFit
{


  // a hacky way to avoid ambiguous base classes
  template<unsigned short NP>
    struct AbstractFitFunctionNames_prim  : public AbstractFitFunctionParameterContainer<std::string,NP>
  {

    typedef AbstractFitFunctionParameterContainer<std::string,NP> AbsBase;

    AbstractFitFunctionNames_prim(void)
      :AbsBase()
    { }

    AbstractFitFunctionNames_prim(const std::string &s)
      :AbsBase(s)
    { }

    ~AbstractFitFunctionNames_prim(void) {}

    virtual void setParNames_prim(const std::vector<std::string> &names)
    {
      AbsBase::set(names); 
    }

    virtual void setParName_prim(const int parNum , const std::string &name)
    {
      AbsBase::set(parNum,name); 
    }

    virtual std::vector<std::string> getParNames_prim(void) const {return AbsBase::get();}

    virtual std::string getParName_prim(const int parNum) const
    {
      return AbsBase::get(parNum); 
    }

    virtual int getParNum_prim(const std::string &name) const 
    {
      return AbsBase::match(name); 
    }

  };


  // relabel names from above


  template<unsigned short NP>
    struct AbstractFitFunctionNames 
    : public AbstractFitFunctionNames_prim<NP>
    {

      typedef AbstractFitFunctionNames_prim<NP> AbsBase;


      AbstractFitFunctionNames(void)
        :AbsBase()
      { }

      AbstractFitFunctionNames(const std::string &s)
        :AbsBase(s)
      { }

      ~AbstractFitFunctionNames(void) {}

      virtual void setParNames(const std::vector<std::string> &names)
      {
        AbsBase::setParNames_prim(names); 
      }

      virtual void setParName(const int parNum , const std::string &name)
      {
        AbsBase::setParName_prim(parNum,name); 
      }

      virtual std::vector<std::string> getParNames(void) const 
      {
        return AbsBase::getParNames_prim();
      }

      virtual std::string getParName(const int parNum) const
      {
        return AbsBase::getParName_prim(parNum); 
      }

      virtual int getParNum(const std::string &name) const 
      {
        return AbsBase::getParNum_prim(name); 
      }

    };


} // jFit











#endif /* ABSTRACT_FUNCTION_PARAMETER_NAMES_H */
