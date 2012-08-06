/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : axis_plotter.cc

 * Purpose :

 * Creation Date : 31-07-2012

 * Last Modified : Wed Aug  1 16:17:23 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "axis_plotter.h"
#include "ensem/ensem.h"
#include "jackFitter/plot.h"
#include "semble/semble_file_management.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/utils/splash.h"
#include <sstream>



namespace radmat
{




  template<>
    void MakeAxisPlots::plot(const ThreePointCorrelator<double> &C3pt,
        const std::string &path, 
        const std::string &fname) const

    {

      SEMBLE::SEMBLEIO::makeDirectoryPath(path);

      ENSEM::EnsemVectorReal cor = C3pt.C3pt;
      std::stringstream axis, jack;

      axis << fname << ".ax";
      jack << fname << ".jack";


      ENSEM::write(jack.str(),cor);

      AxisPlot plot;
      plot.addEnsemData(cor,"\\sq",1);

      plot.sendToFile(axis.str());

    }


  template<>
    void MakeAxisPlots::plot(const ThreePointCorrelator<std::complex<double> > &C3pt,
        const std::string &path,
        const std::string &fname) const
    {
      SEMBLE::SEMBLEIO::makeDirectoryPath(path);

      ENSEM::EnsemVectorComplex tmp;
      tmp = C3pt.C3pt;


      std::stringstream ss_rawj , ss_rawx, ss_real, ss_cplx;

      ss_rawj << fname << "_raw.jack";
      ss_rawx << fname << "_raw.ax" ;
      ss_real << fname << "_real.ax";
      ss_cplx << fname << "_cplx.ax";

      write(ss_rawj.str(),tmp);

      AxisPlot plot,Rplot,Cplot;
      ENSEM::EnsemVectorReal Creal,CComplex;

      Creal = ENSEM::real(tmp);
      CComplex = ENSEM::imag(tmp);

      Rplot.addEnsemData(Creal,"\\sq",1);
      Cplot.addEnsemData(CComplex,"\\sq",2);
      plot.addEnsemData(Creal,"\\sq",1);
      plot.addEnsemData(CComplex,"\\sq",2);

      Rplot.sendToFile(std::string(ss_real.str()));
      Cplot.sendToFile(std::string(ss_cplx.str()));
      plot.sendToFile(std::string(ss_rawx.str()));


    }



} // namespace radmat
