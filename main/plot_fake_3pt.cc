/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : plot_fake_3pt.cc

 * Purpose :

 * Creation Date : 16-07-2012

 * Last Modified : Wed Jul 18 16:22:54 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/fake_data/fake_3pt_function.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "semble/semble_meta.h"
#include <complex>

#include "radmat/utils/splash.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"

#include "jackFitter/plot.h"


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



  ADAT::Handle<FakeDataInputs_p<std::complex<double> > > 
    orig = generateOriginalInputs<std::complex<double> >(fakeinikeys);
  ADAT::Handle<FakeDataInputs<std::complex<double> > > input = copyFakeInput(orig);
  applyZSuppression(input->working);
  applyDispersion(input->working,fakeinikeys.dataProps.momenta[0]);

  Fake3ptCorr<std::complex<double> > foo(input,0,0,0,fakeinikeys.dataProps.momenta[0]);

  std::string raw = std::string("fake_3pt_raw.jack");
  write(raw,foo.get3pt());

  AxisPlot plot,Rplot,Cplot;
  ENSEM::EnsemVectorReal Creal,CComplex;

  Creal = ENSEM::real(foo.get3pt());
  CComplex = ENSEM::imag(foo.get3pt());

  Rplot.addEnsemData(Creal,"\\sq",1);
  Cplot.addEnsemData(CComplex,"\\sq",2);
  plot.addEnsemData(Creal,"\\sq",1);
  plot.addEnsemData(CComplex,"\\sq",2);

  Rplot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);
  Cplot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);
  plot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);


  Rplot.sendToFile(std::string("fake_3pt_real.ax"));
  Cplot.sendToFile(std::string("fake_3pt_imag.ax"));
  plot.sendToFile(std::string("fake_3pt_raw.ax"));


  return 0;
}











