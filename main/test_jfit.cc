/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : test_jfit.cc

 * Purpose :

 * Creation Date : 07-10-2013

 * Last Modified : Mon 14 Oct 2013 05:29:04 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "radmat/jfit/abstract_function_parameter_names.h"
#include "radmat/jfit/abstract_function_default_values.h"
#include "radmat/jfit/abstract_function.h"
#include <iostream>



using namespace jFit; 

  int 
main (void) 
{

  AbstractFitFunctionNames<64> foo;
  int num = 12; 
  std::string n = "I am a foo"; 
  foo.setParName(num,n); 
  std::cout << foo.getParName(num) << " " << num << " " << foo.getParNum(n) << std::endl; 
  

  AbstractFitFunctionDefaultValues<double,64> foo2;
  int num2 = 0; 
  std::string n2 = "I am a foo2" ;
  double v2 = 3.1415;
  foo2.setParName(num2,n2); 
  foo2.setDefaultParValue(n2, v2); 
  std::cout << foo2.getDefaultParValue("I am a foo2") << std::endl;  

  AbstractFitFunction<double,double,64,double,12> func;
  double t0 = 12.21; 
  std::string t = "I am a func"; 
  int pnum = 27;

  func.setParName(pnum,t); 
  std::cout << func.getParNum(t) << std::endl; 
  func.setParamLowerLimit(t,t0); 
  if ( ! func.isParamLowerLimited(pnum) ) 
  {
    std::cerr << "bug up in her" << std::endl; 
    exit (99); 
  } 
  
  std::cout <<func.getParamLowerLimit(t) << std::endl; 

  return 0; 
}
