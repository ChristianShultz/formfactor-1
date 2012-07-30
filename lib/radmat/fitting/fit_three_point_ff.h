#ifndef FIT_THREE_POINT_FF_H_H_GUARD
#define FIT_THREE_POINT_FF_H_H_GUARD


#include "jackFitter/jackknife_fitter.h"
#include "jackFitter/ensem_data.h"
#include "adat/handle.h"
#include "ensem/ensem.h"
#include <string>



namespace radmat
{


  struct Fit3PtFF_Base
  {

    Fit3PtFF_Base(const EnsemData &data, 
        const ADAT::Handle<FitComparator> &fitComp,
        const double noiseRatioCutoff,
        const int minTSlices=4);

    virtual ~Fit3PtFF_Base(void) {}  // do nothing but should be virtual 


    // methods that are common:q

    virtual void saveFitPlot(const std::string &filename) const;
    virtual std::string getFitPlotString(void) const;
    virtual std::string getFitSummary(void) const;
    virtual ENSEM::EnsemReal getFF(void) const;
    virtual double getChisq(void) const;
    virtual int getNDoF(void) const;
    virtual std::string getFitName(void) const;

    // interesting method
    protected:
    virtual void doFit(void) = 0;
    virtual void setFF(const ENSEM::EnsemReal &FF);
    virtual void saveFitInternal(const double xmin,const double xmax);

    ADAT::Handle<FitComparator> m_fitComp;
    double m_noiseRatioCutoff;
    int m_minTSlices;

    std::string m_fit_summary;
    std::string m_axis_plot;
    double m_chisq;
    int m_nDoF;
    std::string m_best_fit_name;

    JackFitLog m_fits;

    ENSEM::EnsemReal m_FF;

    private:
    Fit3PtFF_Base(void);

  };




  struct Fit3PtFF_Constant : public Fit3PtFF_Base
  {
    Fit3PtFF_Constant(const EnsemData &data, 
        const ADAT::Handle<FitComparator> &fitComp, // hate calling these const
        const double noiseRatioCutoff, 
        const int minTSlices=4);

    // interesting method
    protected:
    virtual void doFit(void); 
  };












}


#endif
