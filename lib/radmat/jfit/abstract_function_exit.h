#ifndef ABSTRACT_FUNCTION_EXIT_H
#define ABSTRACT_FUNCTION_EXIT_H

#include <iostream>
#include <stdlib.h>


namespace jFit
{
  struct AbstractFitFunctionExit
  {
    virtual void exit_message(const char *prettyFunction, 
        const char *file, 
        const int line, 
        const char *msg) const
    {
      std::cerr << "Error: " << prettyFunction << " at " 
        << file << ":" << line
        << "\n" << msg << std::endl; 
      exit(EXIT_FAILURE);   
    } 

    virtual void exit_message(const char *prettyFunction, 
        const char *file, 
        const int line, 
        const std::string &msg) const
    { 
      exit_message(prettyFunction,file,line,msg.c_str()); 
    }

    virtual void boolean_exit_message(const bool bin, 
        const char *prettyFunction, 
        const char *file, 
        const int line, 
        const char *msg) const
    {
      if ( bin ) 
        exit_message(prettyFunction,file,line,msg); 
    }

    virtual void boolean_exit_message(const bool bin, 
        const char *prettyFunction, 
        const char *file, 
        const int line, 
        const std::string &msg) const
    {
      if ( bin )
        exit_message(prettyFunction,file,line,msg); 
    }

  };
} // jFit


#endif /* ABSTRACT_FUNCTION_EXIT_H */
