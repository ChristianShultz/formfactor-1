#ifndef ABSTRACT_FUNCTION_H
#define ABSTRACT_FUNCTION_H

#include <vector>
#include <string>
#include "abstract_function_parameter_names.h"
#include "abstract_function_default_errors.h"
#include "abstract_function_default_values.h"
#include "abstract_function_parameter_fixing.h"
#include "abstract_function_parameter_limits.h"


namespace jFit
{

  template<typename RType,      // return type
    typename PType,             // parameter type <vector>
    unsigned short NP,          // npars
    typename IType,             // input type <vector>
    unsigned short NI >         // ninputs
      struct AbstractFitFunction 
      :public AbstractFitFunctionNames<NP>,
      public AbstractFitFunctionDefaultErrors<PType,NP>, 
      public AbstractFitFunctionDefaultValues<PType,NP>, 
      public AbstractFitFunctionParameterFixing<NP>, 
      public AbstractFitFunctionParameterLimits<PType,NP> 
  {

    AbstractFitFunction(void)
      : AbstractFitFunctionNames<NP>::AbstractFitFunctionNames(),
      AbstractFitFunctionDefaultErrors<PType,NP>::AbstractFitFunctionDefaultErrors(), 
      AbstractFitFunctionDefaultValues<PType,NP>::AbstractFitFunctionDefaultValues(), 
      AbstractFitFunctionParameterFixing<NP>::AbstractFitFunctionParameterFixing(), 
      AbstractFitFunctionParameterLimits<PType,NP>::AbstractFitFunctionParameterLimits() 
    {}
    ~AbstractFitFunction(void) {};


    virtual RType operator()(const std::vector<IType> &) const {return RType(1);};
    virtual std::string getFitType(void) const {return "default";}; 
    /*    virtual RType operator()(const std::vector<IType> &) const = 0;
          virtual std::string getFitType(void) const = 0; */

    virtual unsigned short getNPars(void) const {return NP;}

    // override to avoid ambiguious base
    virtual void setParNames(const std::vector<std::string> &names)
    {
      AbstractFitFunctionNames<NP>::setParNames(names); 
      AbstractFitFunctionDefaultErrors<PType,NP>::setParNames(names);  
      AbstractFitFunctionDefaultValues<PType,NP>::setParNames(names); 
      AbstractFitFunctionParameterFixing<NP>::setParNames(names); 
      AbstractFitFunctionParameterLimits<PType,NP>::setParNames(names); 
    }

    virtual void setParName(const int pnum, const std::string &n)
    {
      AbstractFitFunctionNames<NP>::setParName(pnum,n); 
      AbstractFitFunctionDefaultErrors<PType,NP>::setParName(pnum,n);  
      AbstractFitFunctionDefaultValues<PType,NP>::setParName(pnum,n); 
      AbstractFitFunctionParameterFixing<NP>::setParName(pnum,n); 
      AbstractFitFunctionParameterLimits<PType,NP>::setParName(pnum,n); 
    }

  };



} // jFit

#endif /* ABSTRACT_FUNCTION_H */
