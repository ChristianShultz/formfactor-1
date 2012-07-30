// fit_three_point_ff.cc -
//
// Friday, July 27 2012
//

#include"fit_three_point_ff.h"
#include <vector>
#include "radmat/utils/pow2assert.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include "adat/handle.h"
#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/jackknife_fitter.h"
#include "radmat_fit_forms.h"
#include "radmat/utils/splash.h"

namespace radmat
{

  ////////////////////////////////////////
  //       **  Fit3PtFF_Base  **        //
  ////////////////////////////////////////

  Fit3PtFF_Base::Fit3PtFF_Base(const EnsemData &data, 
      const ADAT::Handle<FitComparator> &fitComp,
      const double noiseRatioCutoff,
      const int minTSlices /* =4 */)
    : m_fitComp(fitComp) ,
    m_noiseRatioCutoff(noiseRatioCutoff) ,
    m_minTSlices(minTSlices),
    m_fits(data)
  {}

  void Fit3PtFF_Base::saveFitPlot(const std::string &file) const
  {
    std::ofstream out;
    out.open(file.c_str());
    out << m_axis_plot;
    out.close();
  }

  std::string Fit3PtFF_Base::getFitPlotString(void) const
  {
    return m_axis_plot;
  }

  std::string Fit3PtFF_Base::getFitSummary(void) const
  {
    return m_fit_summary;
  }

  ENSEM::EnsemReal  Fit3PtFF_Base::getFF(void) const
  {
    return m_FF;
  }

  double  Fit3PtFF_Base::getChisq(void) const
  {
    return m_chisq;
  }

  int  Fit3PtFF_Base::getNDoF(void) const
  {
    return m_nDoF;
  }

  std::string  Fit3PtFF_Base::getFitName(void) const
  {
    return m_best_fit_name;
  }

  void Fit3PtFF_Base::setFF(const ENSEM::EnsemReal &FF) 
  {
    m_FF = FF;
  }

  void Fit3PtFF_Base::saveFitInternal(const double xmin, const double xmax)
  {
    int rank; // this isn't actually used in the next method.. ?? 
    FitDescriptor best_fit_descriptor = m_fits.getBestJackFit(*m_fitComp,rank);
    if(best_fit_descriptor.fitname == "FAILED")
    {
      SPLASH("fit failed");

      ENSEM::EnsemReal dum;
      dum.resize(m_fits.getEnsemData().getNBins());
      dum = ENSEM::Real(0.);
      setFF(dum);
      m_chisq = 1.0e10;
      m_nDoF = 1;
      m_best_fit_name = std::string("FAILED");
      m_fit_summary = std::string("FAILED -  set Z to 0.0");
      //no plot !!!
    }
    else
    {
      setFF((m_fits.getFit(best_fit_descriptor)).getJackFitParValue("FF"));
      m_chisq = (m_fits.getFit(best_fit_descriptor)).getJackChisq();
      m_nDoF = (m_fits.getFit(best_fit_descriptor)).getNDoF();
      m_best_fit_name = best_fit_descriptor.fitname;

      // axis plot
      std::stringstream label;
      label <<  "\\gx\\sp2\\ep/N\\sbdof\\eb=" << std::setprecision(2) << m_chisq<< "/" << m_nDoF
        << "; FF=" << fixed << std::setprecision(4) << ENSEM::toDouble(ENSEM::mean(m_FF)) << "\\+-" 
        <<  std::setprecision(4) << ENSEM::toDouble(ENSEM::sqrt(ENSEM::variance(m_FF)));


      m_axis_plot = (m_fits.getFit(best_fit_descriptor)).makeJackFitPlotAxis(xmin,xmax,label.str());

    }


  }



  ////////////////////////////////////////
  //      **  Fit3PtFF_Const  **        //
  ////////////////////////////////////////


  Fit3PtFF_Constant::Fit3PtFF_Constant(const EnsemData &data, 
      const ADAT::Handle<FitComparator> &fitComp,
      const double noiseRatioCutoff,
      const int minTSlices /* =4 */)
    : Fit3PtFF_Base(data,fitComp,noiseRatioCutoff,minTSlices)
  {
    doFit();
  }

  void Fit3PtFF_Constant::doFit(void)
  {
    // hide the baddies
    m_fits.getEnsemData().hideDataAboveYErrRat(m_noiseRatioCutoff);

    std::vector<double> all_time = m_fits.getEnsemData().getAllXData();

    // check that there are enough time slices to fit  
    POW2_ASSERT((int)all_time.size() >= m_minTSlices);

    // find the ranges
    int t_high = int( *std::max_element(all_time.begin(),all_time.end()) );
    int t_low = int(*std::min_element(all_time.begin(),all_time.end()) );

    ADAT::Handle<FitFunction> m_Constant(new radmat::Constant());
    POW2_ASSERT(&*m_Constant);

    ENSEM::EnsemReal m_default_par;
    double tmid = double(t_low-t_high)/2.;
    m_default_par = m_fits.getEnsemData().getYUsingNearestX(tmid);

    m_Constant->setDefaultParValue("FF",ENSEM::toDouble(ENSEM::mean(m_default_par)));
    m_Constant->setDefaultParError("FF",ENSEM::toDouble(ENSEM::sqrt(ENSEM::variance(m_default_par))));

    // vary start and end points searching for optimal fit
    for(int t_small = t_low; t_high - t_small >= m_minTSlices;  ++t_small)
      for(int t_big = t_high; t_big - t_small >= m_minTSlices; --t_big)
      {
        m_fits.getEnsemData().showAll();
        m_fits.getEnsemData().hideDataAboveX(t_big + 0.1);
        m_fits.getEnsemData().hideDataBelowX(t_small - 0.1);
        std::stringstream ss; 
        ss << "const tmin = " << t_small << " tmax = " << t_big;

        // could be baddies in the range
        if(m_fits.getEnsemData().getNData() >= m_minTSlices)
          m_fits.addFit(ss.str(),m_Constant);

      }

    saveFitInternal(double(t_low) -0.5, double(t_high) + 0.5);

    m_fit_summary = std::string("CHRISTIAN NEEDS TO IMPLEMENT FIT SUMMARIES STILL");

  }








}
