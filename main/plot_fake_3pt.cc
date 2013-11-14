/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : plot_fake_3pt.cc

 * Purpose :

 * Creation Date : 16-07-2012

 * Last Modified : Fri Nov  9 11:23:51 2012

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

#include "radmat/load_data/load_fake_data.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/load_data/three_point.h"
#include "radmat/fake_data/fake_data_ini.h"
#include "radmat/llsq/llsq_driver.h"
#include "radmat/llsq/llsq_q2_pack.h"
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


  /*
     rHandle<FakeDataInputs_p<std::complex<double> > > 
     orig = generateOriginalInputs<std::complex<double> >(fakeinikeys);
     rHandle<FakeDataInputs<std::complex<double> > > input = copyFakeInput(orig);
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

   */



  LoadFake3pt<std::complex<double> > loader(fakeinikeys);

  std::vector<ThreePointCorrelator<std::complex<double> > > C3 = loader.genData();

  int ct(0);
  std::vector<ThreePointCorrelator<std::complex<double> > >::const_iterator it;
  for(it = C3.begin(); it != C3.end(); it++,ct++)
  {

    ENSEM::EnsemVectorComplex tmp = it->C3pt;

    std::stringstream ss, ss_rawj , ss_rawx, ss_real, ss_cplx;
    ss << "fake_3pt_" << ct;

    ss_rawj << ss.str() << "_raw.jack";
    ss_rawx << ss.str() << "_raw.ax" ;
    ss_real << ss.str() << "_real.ax";
    ss_cplx << ss.str() << "_cplx.ax";

    std::string raw = std::string(ss_rawj.str());
    write(raw,tmp);

    AxisPlot plot,Rplot,Cplot;
    ENSEM::EnsemVectorReal Creal,CComplex;

    Creal = ENSEM::real(tmp);
    CComplex = ENSEM::imag(tmp);

    Rplot.addEnsemData(Creal,"\\sq",1);
    Cplot.addEnsemData(CComplex,"\\sq",2);
    plot.addEnsemData(Creal,"\\sq",1);
    plot.addEnsemData(CComplex,"\\sq",2);

    Rplot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);
    Cplot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);
    plot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);


    Rplot.sendToFile(std::string(ss_real.str()));
    Cplot.sendToFile(std::string(ss_cplx.str()));
    plot.sendToFile(std::string(ss_rawx.str()));
  }



  BuildQ2Packs<std::complex<double> > Q2Builder;
  Q2Builder.load(C3);
  Q2Builder.strip_propagation_factor();

  std::vector<rHandle<LLSQDataPointQ2Pack> > q2_packs = Q2Builder.getQ2Packs();
  std::vector<rHandle<LLSQDataPointQ2Pack> >::const_iterator q2_pack_iterator;

  int major_index(0);
  for(q2_pack_iterator = q2_packs.begin(); q2_pack_iterator != q2_packs.end(); ++q2_pack_iterator, ++major_index)
  {
    LLSQRet_ff_Q2Pack<std::complex<double> > foobar;
    LLSQDriver_t<std::complex<double> > llsq_driver(std::string("SVDMakeSquare"));
    foobar = *llsq_driver(*q2_pack_iterator);

    LLSQRet_ff_Q2Pack<std::complex<double> >::const_iterator it;


    for(it = foobar.begin(); it != foobar.end(); ++it)
    {

      ENSEM::EnsemVectorComplex tmp = it->second;
      std::stringstream ss, ss_rawj, ss_rawx, ss_cplx, ss_real;
      ss << "FF_Q2_" << major_index << "_ffnum_" << it->first;

      ss_rawj << ss.str() << "_raw.jack";
      ss_rawx << ss.str() << "_raw.ax" ;
      ss_real << ss.str() << "_real.ax";
      ss_cplx << ss.str() << "_cplx.ax";


      std::string raw = std::string(ss_rawj.str());
      write(raw,tmp);

      AxisPlot plot,Rplot,Cplot;
      ENSEM::EnsemVectorReal Creal,CComplex;

      Creal = ENSEM::real(tmp);
      CComplex = ENSEM::imag(tmp);

      Rplot.addEnsemData(Creal,"\\sq",1);
      Cplot.addEnsemData(CComplex,"\\sq",2);
      plot.addEnsemData(Creal,"\\sq",1);
      plot.addEnsemData(CComplex,"\\sq",2);

      Rplot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);
      Cplot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);
      plot.setXRange(fakeinikeys.timeProps.tsource -2, fakeinikeys.timeProps.tsink +2);


      Rplot.sendToFile(std::string(ss_real.str()));
      Cplot.sendToFile(std::string(ss_cplx.str()));
      plot.sendToFile(std::string(ss_rawx.str()));

    }
  }





  return 0;
}











