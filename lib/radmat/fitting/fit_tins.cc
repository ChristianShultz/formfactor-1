/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : fit_tins.cc

 * Purpose :

 * Creation Date : 01-08-2012

 * Last Modified : Mon 30 Sep 2013 05:11:52 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "fit_tins.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "semble/semble_vector.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "jackFitter/three_point_fit_forms.h"
#include "jackFitter/plot.h"
#include "adat/handle.h"
#include <string>
#include <math.h>
#include <vector>
#include <iostream>
#include <complex>
#include "radmat/utils/splash.h"
#include "radmat/utils/pow2assert.h"




namespace radmat
{


  namespace 
  {
    inline double rad(const double deg)
    {
      return 3.14159 * deg/180.;
    }

    // this is a bit hackey since it includes contact terms and the ones where the insertion has run out past the source/sink
    ENSEM::EnsemVectorReal convertToReal(const TinsFitter &fitter, const ENSEM::EnsemVectorComplex & in, const int tlow, const int thigh)
    {
      ENSEM::EnsemVectorReal real, imag;

      double mean_phase(0); 

      real = ENSEM::real(in);
      imag = ENSEM::imag(in);


      for(int i = tlow; i < thigh; ++i)
      {
        double rl,im;
        rl = SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(real,i)));
        im = SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(imag,i)));
        
        double phase = std::arg(std::complex<double>(rl,im)); 
        
        
        // std::cout << __func__ << " t = " << i << " phase = " << phase << " (" << phase*180./3.14159 <<" deg)"<< std::endl;

      
        if(phase < -3.*3.14159/4.)
          phase = - phase; 
    

        //  std::cout << __func__ << " t = " << i << " phase = " << phase << " (" << phase*180./3.14159 <<" deg)"<< std::endl;

        mean_phase += phase; 
      }


      double phase = mean_phase/double(thigh - tlow); 



      if( (phase < 0.174528) && (phase > -0.174708) ) // +/- 10 degree about 0 in rad
        return real;
      else if( (phase > 1.39622) && (phase < 1.74528)) // 90deg
        return imag;
      else if( (phase > 2.96697) || (phase < -2.96715))  //180deg // this is atan2 specific, it returns (-pi,pi)
        return real;
      else if( (phase > -1.74546) && (phase < -1.3964)) // 270 de/doFit
        return imag;
      else
      {
        std::cout << "The calculated phase was " << phase*180./3.14159 << " (deg)" << std::endl;
        std::cout << "for Q2 = " << SEMBLE::toScalar(ENSEM::mean(fitter.getQ2())) << std::endl;  
        SPLASH("check bad_corr.jack, bad_corr.ax for the correlator"); 

        AxisPlot plot; 
        plot.addEnsemData(ENSEM::real(in),"\\sq",1);
        plot.addEnsemData(ENSEM::imag(in),"\\sq",2);
        plot.sendToFile("bad_corr.ax");

        ENSEM::write("bad_corr.jack",in); 

        SPLASH("An error occured while trying to pull the real/imag part of the solution vector,exiting.");
        exit(1);
      }
    }


  } // anonymous


  template<>
    void TinsFitter::fit<double>(const std::string &fname, 
        const LLSQRet_ff_Q2Pack<double> &pack,
        const ThreePointComparatorProps_t &fitProps)
    {
      const unsigned int sz = pack.size();
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,sz);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<double>::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,it->second,it->first,fitProps);

      didFit = true;
    }



  template<>
    void TinsFitter::fit<std::complex<double> >(const std::string &fname , 
        const LLSQRet_ff_Q2Pack<std::complex<double> > &pack,
        const ThreePointComparatorProps_t &fitProps)
    {
      const unsigned int sz = pack.size();
      const int nbins = pack.Q2().size();
      ff.reDim(nbins,sz);
      Q2 = pack.Q2();

      LLSQRet_ff_Q2Pack<std::complex<double> >::const_iterator it;

      for(it = pack.begin(); it != pack.end(); ++it)
        doFit(fname,convertToReal(*this,it->second,fitProps.tlow,fitProps.thigh),it->first,fitProps);

      didFit = true;
    }


  void TinsFitter::doFit(const std::string &filenameBase, 
      const ENSEM::EnsemVectorReal &data, 
      const int ffnum,
      const ThreePointComparatorProps_t &fitProps)
  {
    std::stringstream ss,jack,ax,fit;
    ss << filenameBase;
    fit << ss.str() << "FF_" << ffnum << "_fit.jack";
    jack << ss.str() << "FF_" << ffnum << ".jack";
    ax << ss.str() << "FF_" << ffnum << ".ax";

    std::vector<double> time; 
    const int Lt = data.numElem();

    for(int t = 0; t < Lt; t++)
      time.push_back(t);

    EnsemData corrData(time,data);


    ADAT::Handle<FitComparator> fitComp = constructThreePointFitComparator(fitProps);
    ADAT::Handle<FitThreePoint> fitCorr (new FitThreePoint(corrData,Lt,0,fitComp,fitProps.minTSlice,fitProps.fit_type));

    // send to files
    fitCorr->saveFitPlot(ax.str());
    write(jack.str(),data);

    ff.loadEnsemElement(ffnum,fitCorr->getFF());
    write(fit.str(),fitCorr->getFF()); 

    fitters.insert(std::map<int,ADAT::Handle<FitThreePoint> >::value_type(ffnum,fitCorr));
  }


  void TinsFitter::writeFitLogs(const std::string &path) const
  {
    std::map<int,ADAT::Handle<FitThreePoint> >::const_iterator it;
    for(it = fitters.begin(); it != fitters.end(); ++it)
    {
      std::stringstream ss;
      ss << path << "_F_" << it->first << "_fit_log.log";
      std::ofstream out(ss.str().c_str());
      out << it->second->getFitSummary();
      out.close();
    }

  }

  void TinsFitter::writeFitPlotsWithComponents(const std::string & path) const
  {
    std::map<int,ADAT::Handle<FitThreePoint> >::const_iterator it;
    for(it = fitters.begin(); it != fitters.end(); ++it)
    {
      std::stringstream ss; 
      ss << path <<  "_F_" << it->first << "_component_fit.ax";
      std::ofstream out(ss.str().c_str());
      out << it->second->getFitPlotStringWithComponents(); 
      out.close();
    }

  }


  std::pair<ENSEM::EnsemReal,SEMBLE::SembleVector<double> > TinsFitter::fetchFF(void) const
  {
    POW2_ASSERT(didFit);

    return std::pair<ENSEM::EnsemReal,SEMBLE::SembleVector<double> >(Q2,ff);
  }


  ENSEM::EnsemReal TinsFitter::getFF(const int index) const
  {
    POW2_ASSERT(index < ff.getN());
    return ff.getEnsemElement(index);
  }

  ADAT::Handle<FitThreePoint> TinsFitter::getFit(const int ffnum) const
  {
    POW2_ASSERT(fitters.find(ffnum) != fitters.end());
    return ADAT::Handle<FitThreePoint>(fitters.find(ffnum)->second);
  }


} // namespace radmat
