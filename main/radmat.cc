/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : radmat.cc

* Purpose :

* Creation Date : 10-12-2012

* Last Modified : Mon Dec 10 09:26:33 2012

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "semble/semble_meta.h"
#include <complex>

#include "radmat/utils/splash.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"

#include "jackFitter/plot.h"

#include "radmat/load_data/load_fake_data.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/load_data/three_point.h"
#include "radmat/load_data/build_correlators.h"
#include "radmat/llsq/llsq_driver.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "jackFitter/plot.h"

#include "radmat/driver/radmat_driver.h"
#include "radmat/driver/radmat_driver_props.h"

#include "semble/semble_file_management.h"


using namespace radmat;
using namespace ENSEM;
using namespace ADAT;
using namespace ADATIO;

  int 
main(int argc, char *argv[])
{
  if(argc != 3)
  {
    SPLASH("usage: test_fake_ini : <xmlinifile> <driverProps> ");
    exit(1);
  }

  std::string xmlinifile;
  std::istringstream val(argv[1]);
  val >> xmlinifile;
  std::cout << "Loading xmlinifile: " << xmlinifile << std::endl;

ThreePointCorrIni_t inikeys; 
  try
  {
    XMLReader xml(xmlinifile);
    read(xml,"/inikeys",inikeys);
  }
  catch(std::exception &e)
  {
    std::cout << "std exception: " << e.what();
  }
  catch(std::string &e)
  {
    std::cout << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }
  catch(...)
  {
    SPLASH("An error occured while loading the fake data inifile");
    exit(1);
  }


  std::string xmlinifile2;
  std::istringstream val2(argv[2]);
  val2 >> xmlinifile2;
  std::cout << "Loading xmlinifile: " << xmlinifile2 << std::endl;

  RDriverProps_t driverProps;

  try
  {
    XMLReader xml(xmlinifile2);
    read(xml,"/DriverProps",driverProps);
  }
  catch(std::exception &e)
  {
    std::cout << "std exception: " << e.what();
  }
  catch(std::string &e)
  {
    std::cout << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }
  catch(...)
  {
    SPLASH("An error occured while loading the fake data inifile");
    exit(1);
  }


  std::cout << "Loading correlators and inverting subduction.. " << std::endl; 

  BuildCorrelators corr_builder(inikeys); 



  RDriver<std::complex<double> > driver(driverProps);
  driver.load(corr_builder.build_correlators());
  std::cout << "Solving llsq.." << std::endl;

  driver.solve_llsq();


  for(int i = 0; i < driver.getNumFFs(); ++i)
  {
    std::cout << "Solving F_" << i << "(Q2).." << std::endl;

    AxisPlot plot = driver.getPlot(i);

    std::stringstream filename;
    filename << SEMBLE::SEMBLEIO::getPath() << "FormFactors/";

    SEMBLE::SEMBLEIO::makeDirectoryPath(filename.str());

    filename << "FF_" << i << "_of_Q2.ax";

    plot.sendToFile(filename.str());
  }


  return 0;
}
