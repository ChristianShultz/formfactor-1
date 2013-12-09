/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : radmat_ff.cc

* Purpose :

* Creation Date : 06-12-2013

* Last Modified : Fri 06 Dec 2013 08:55:10 AM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include <string>
#include <sstream>
#include <iostream>
#include <utility>
#include "radmat/driver/radmat_driver.h"
std::string 
usage(int argc, char *argv[])
{
  if(argc != 2)
  {
    std::cerr << "usage: radmat_ff : <ini> " << std::endl;
    exit(1); 
  }

  std::string xmlini;
  std::istringstream val(argv[1]);
  val >> xmlini;


  return xmlini;
}


