// fake_radmat.cc -
//
// Wednesday, August  1 2012
//



#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/fake_data/fake_3pt_function.h"
#include "radmat/fake_data/fake_data_ini.h"
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
#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/llsq/llsq_driver.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "jackFitter/plot.h"

#include "radmat/driver/radmat_driver.h"

#include "semble/semble_file_management.h"


using namespace radmat;
using namespace ENSEM;
using namespace ADAT;
using namespace ADATIO;

  int 
main(int argc, char *argv[])
{
  if(argc != 2)
  {
    SPLASH("usage: test_fake_ini : <xmlinifile> ");
    exit(1);
  }

  std::string xmlinifile;
  std::istringstream val(argv[1]);
  val >> xmlinifile;
  std::cout << "Loading xmlinifile: " << xmlinifile << std::endl;

  FakeDataIni_t fakeinikeys;

  try
  {
    XMLReader xml(xmlinifile);
    read(xml,"/FakeDataIni",fakeinikeys);
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


  LoadFake3pt<std::complex<double> > loader(fakeinikeys);

    
  std::vector<ThreePointCorrelator<std::complex<double> > > C3 = loader.genData();

  std::vector<LoadFake3pt<std::complex<double> >::ffFunction > formfactors;
  formfactors = loader.get_FF_inputs();

  RDriver<std::complex<double> > driver;
  driver.load(C3);
  driver.solve_llsq();

  
  for(int i = 0; i < driver.getNumFFs(); ++i)
    {
      AxisPlot plot = driver.getPlot(i);
      std::pair<double,double> range = driver.getQ2Range();
      std::vector<double> q2, ff;

      range.first = range.first - 0.1;
      range.second = range.second +0.1;
  

      POW2_ASSERT(range.first <= range.second); // loop sanity

      for(double dq2 = range.first; dq2 < range.second; dq2 += (range.second - range.first)/double(50))
        {
          q2.push_back(dq2);
          ff.push_back(formfactors[i](dq2));
        }

      plot.addLineData(q2,ff,2);

      std::stringstream filename;
      filename << SEMBLE::SEMBLEIO::getPath() << "FormFactors/";
      
      SEMBLE::SEMBLEIO::makeDirectoryPath(filename.str());

      filename << "FF_" << i << "_of_Q2.ax";
        
      plot.sendToFile(filename.str());
    }


  return 0;
}
