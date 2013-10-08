#ifndef ABSTRACT_FUNCTION_PARAMETER_NAMES_H
#define ABSTRACT_FUNCTION_PARAMETER_NAMES_H 


#include <string>
#include <vector>
#include <sstream>
#include "abstract_function_exit.h"


namespace jFit
{

  template<unsigned short NP>
    struct AbstractFitFunctionNames  : public AbstractFitFunctionExit
  {
    AbstractFitFunctionNames(void)
    {
      parNames.resize(NP,"default"); 
    }

    ~AbstractFitFunctionNames(void) {}

    virtual void setParNames(const std::vector<std::string> &names)
    {
      parNames = names; 
    }

    virtual void setParName(const int parNum , const std::string &name)
    {
      boolean_exit_message( parNum > NP ,
          __PRETTY_FUNCTION__, 
          __FILE__ , 
          __LINE__ , 
          "parNum was larger than nPar");

      parNames[parNum] = name;  
    }

    virtual std::vector<std::string> getParNames(void) const {return parNames;}

    virtual std::string getParName(const int parNum) const
    {
      boolean_exit_message( parNum > NP ,
          __PRETTY_FUNCTION__, 
          __FILE__ , 
          __LINE__ , 
          "parNum was larger than nPar");

      return parNames[parNum]; 
    }

    virtual int getParNum(const std::string &name) const 
    {
      for (int i = 0; i < NP; ++i)
        if ( parNames[i] == name ) 
          return i;

      std::stringstream ss; 
        ss << "unable to match " << name << std::endl ; 

      exit_message(__PRETTY_FUNCTION__, 
          __FILE__ , 
          __LINE__ , 
          ss.str() );

    }

    std::vector<std::string> parNames; 
  };


} // jFit











#endif /* ABSTRACT_FUNCTION_PARAMETER_NAMES_H */
