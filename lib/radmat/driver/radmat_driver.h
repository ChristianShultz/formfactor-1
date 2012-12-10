#ifndef RADMAT_DRIVER_H_H_GUARD
#define RADMAT_DRIVER_H_H_GUARD

#include "radmat/load_data/three_point.h"
#include "radmat/load_data/build_q2_packs.h"
#include "radmat/load_data/build_correlators.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "radmat/fitting/axis_plotter.h"
#include "radmat/fitting/fit_tins.h"
#include "radmat_driver_props.h"
#include "adat/handle.h"
#include "semble/semble_file_management.h"
#include "semble/semble_meta.h"
#include "jackFitter/plot.h"
#include <vector>
#include <string>
#include <sstream>
#include <string>
#include <algorithm>
#include <utility>
#include <iostream>

#define USE_OMP_RADAMAT_DRIVER_SOLVE_LLSQ

#ifdef  USE_OMP_RADAMAT_DRIVER_SOLVE_LLSQ
#include <omp.h>
#endif


namespace radmat
{

  template<typename T> 
    struct RDriver
    {
      RDriver(void);
      RDriver(const RDriverProps_t &driverProps); 
      RDriver(const RDriver<T> &o);
      ~RDriver(void){clear();}
      RDriver<T>& operator=(const RDriver<T> &o);


      void load(const std::vector<ThreePointCorrelator<T> > &corrs);
      void load(const std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > &q2_packs);

      void solve_llsq(void);

      int getNumFFs(void) const;

      AxisPlot getPlot(const int ffnum) const;

      std::pair<double,double> getQ2Range(void) const;

      private:

      void clear(void);
      void makeQ2Plots(void);

      RDriverProps_t m_driverProps; 
      

      bool haveQ2Packs;
      std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > m_q2_packs;

      bool solvedForFF;
      std::vector<ADAT::Handle<LLSQRet_ff_Q2Pack<T> > > m_q2_ff_packs;
      std::vector<TinsFitter> m_fits;

      bool madePlots;
      std::vector<AxisPlot> m_ff_of_q2;

    };

#if 0 // -- hide ctor
  template<typename T>
    RDriver<T>::RDriver(void)
    : haveQ2Packs(false), solvedForFF(false), madePlots(false)
    {  }
#endif

  template<typename T>
    RDriver<T>::RDriver(const RDriverProps_t &driverProps)
    : m_driverProps(driverProps), haveQ2Packs(false), solvedForFF(false), madePlots(false)
    {  }


  template<typename T>
    RDriver<T>::RDriver(const RDriver<T> &o)
    : m_driverProps(o.m_driverProps),haveQ2Packs(o.haveQ2Packs) , m_q2_packs(o.m_q2_packs), 
    solvedForFF(o.solvedForFF), m_q2_ff_packs(o.m_q2_ff_packs), m_fits(o.m_fits) , 
    madePlots(o.madePlots), m_ff_of_q2(o.m_ff_of_q2)
  {  }


  template<typename T>
    RDriver<T>& RDriver<T>::operator=(const RDriver<T> &o)
    {
      if(this != &o)
      {
        m_driverProps = o.m_driverProps; 
        haveQ2Packs = o.haveQ2Packs;
        m_q2_packs = o.m_q2_packs;
        solvedForFF = o.solvedForFF;
        m_q2_ff_packs = o.m_q2_ff_packs; 
        m_fits = o.m_fits;
        madePlots = o.madePlots;
        m_ff_of_q2 = o.m_ff_of_q2;
      }
      return *this;
    }


  template<typename T>
    void RDriver<T>::load(const std::vector<ThreePointCorrelator<T> > &corrs) 
    {

      clear();

      typename std::vector<ThreePointCorrelator<T> >::const_iterator it;
      MakeAxisPlots m_plotter;    

      for(it = corrs.begin(); it != corrs.end(); ++it)
      {
        // convert Q2 to a string format
        double Q2 = SEMBLE::toScalar(ENSEM::mean(it->Q2));
        std::stringstream ssQ2;
        ssQ2 << Q2;
        std::string sQ2 = ssQ2.str();
        std::replace(sQ2.begin(),sQ2.end(),'.','p');
        std::replace(sQ2.begin(),sQ2.end(),'-','m');

        std::string path = SEMBLE::SEMBLEIO::getPath();
        path += std::string("ThreePointCorrelators");
        std::stringstream fname;
        fname << path << "/C3pt_l" << it->lorentz  << "_hf" << it->hel_sink 
          << "_hi" << it->hel_source << "_ID_" << it->elemIDBase << "_pf" << it->mom.momSink[0] << "_"
          << it->mom.momSink[1] << "_" << it->mom.momSink[2] << "_pi" << it->mom.momSource[0] << "_"
          << it->mom.momSource[1] << "_" << it->mom.momSource[2] << "_Q2_" << sQ2;

        m_plotter.plot(*it,path,fname.str());
      }


      BuildQ2Packs<T> Q2Builder;
      Q2Builder.load(corrs);
      Q2Builder.strip_propagation_factor();
      load(Q2Builder.getQ2Packs());
    }

  template<typename T>
    void RDriver<T>::load(const std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > &q2_packs)
    {
      clear();
      m_q2_packs = q2_packs;
      haveQ2Packs = true;
    }


  template<typename T>
    void RDriver<T>::solve_llsq(void)
    {
      int sz = m_q2_packs.size();
      int pack_index;
      m_q2_ff_packs.resize(sz);

      m_fits.resize(sz);

#ifdef  USE_OMP_RADAMAT_DRIVER_SOLVE_LLSQ
#pragma omp parallel for shared(pack_index, sz)
#endif

      // solve the llsq and then fit out the t_ins dependence
      for(pack_index = 0; pack_index < sz; ++pack_index)
      {
        LLSQDriver_t<T> m_llsq_driver(std::string("SVDMakeSquare"));
        m_q2_ff_packs[pack_index] = m_llsq_driver(m_q2_packs[pack_index]);

        // write out the axis plot of F_n(Q2,t_ins)
        std::string path = SEMBLE::SEMBLEIO::getPath();
        path += std::string("t_ins_fits");
        SEMBLE::SEMBLEIO::makeDirectoryPath(path);  

        // convert Q2 to a string format
        ENSEM::EnsemReal EQ2 = m_q2_packs[pack_index]->begin()->second.begin()->Q2();
        double Q2 = SEMBLE::toScalar(ENSEM::mean(EQ2));
        std::stringstream ssQ2;
        ssQ2 << Q2;
        std::string sQ2 = ssQ2.str();

        std::replace(sQ2.begin(),sQ2.end(),'.','p');
        std::replace(sQ2.begin(),sQ2.end(),'-','m');     


        // do/write out fits
        std::stringstream ssfname;
        ssfname << path << "/FF_Q2_" << sQ2;
        m_fits[pack_index].fit<T>(ssfname.str(),
            *m_q2_ff_packs[pack_index],
            m_driverProps.threePointComparatorProps);

        // write out fit logs
        std::stringstream log;
        log << path << "/fit_logs";
        SEMBLE::SEMBLEIO::makeDirectoryPath(log.str());
        m_fits[pack_index].writeFitLogs(log.str() + std::string("/FF_Q2_") + sQ2);

        // write out the fits with the components plotted
        std::stringstream componentFits; 
        componentFits << path << "/ComponentFits";
        SEMBLE::SEMBLEIO::makeDirectoryPath(componentFits.str());
        m_fits[pack_index].writeFitPlotsWithComponents(componentFits.str() +  std::string("/FF_Q2_") + sQ2);


      }

      solvedForFF = true;

      makeQ2Plots();
    }


  template<typename T>
    int RDriver<T>::getNumFFs(void) const
    {
      POW2_ASSERT(madePlots);
      return m_ff_of_q2.size();
    }


  template<typename T>
    AxisPlot RDriver<T>::getPlot(const int ff) const
    {
      POW2_ASSERT(ff < getNumFFs());
      return m_ff_of_q2[ff];
    }


  template<typename T>
    void RDriver<T>::makeQ2Plots(void)
    {

      m_ff_of_q2.clear();

      const int nff = m_fits.begin()->nFF();
      const int nQ2s = m_fits.size();
      const int ncfg = m_fits.begin()->getQ2().size();

      std::vector<ENSEM::EnsemReal> Q2s;
      std::vector<double> dQ2s;
      for(int elem = 0; elem < nQ2s; ++elem)
      {
        Q2s.push_back( m_fits[elem].getQ2() );
        dQ2s.push_back(SEMBLE::toScalar(ENSEM::mean(m_fits[elem].getQ2())));
      }

      for(int ff = 0; ff < nff; ++ff)
      {
        ENSEM::EnsemVectorReal FF;
        FF.resize(ncfg);
        FF.resizeObs(nQ2s);

        for(int elem = 0; elem < nQ2s; ++elem)
          ENSEM::pokeObs(FF,m_fits[elem].getFF(ff),elem);

        AxisPlot plot;
        plot.addEnsemData(dQ2s,FF,"\\sq",1);

        // add a label here..

        m_ff_of_q2.push_back(plot);

      }

      madePlots = true;
    }



  template<typename T>
    std::pair<double,double> RDriver<T>::getQ2Range(void) const
    {
      POW2_ASSERT(solvedForFF);
      const int nQ2s = m_fits.size();
      std::vector<double> dQ2s;
      for(int elem = 0; elem < nQ2s; ++elem)
        dQ2s.push_back(SEMBLE::toScalar(ENSEM::mean(m_fits[elem].getQ2())));


      return std::pair<double,double>(*std::min_element(dQ2s.begin(),dQ2s.end()),
          *std::max_element(dQ2s.begin(),dQ2s.end()));
    }



  template<typename T>
    void RDriver<T>::clear(void)
    {
      m_q2_packs.clear();
      haveQ2Packs = false;
      m_q2_ff_packs.clear();
      m_fits.clear();
      solvedForFF = false;
      m_ff_of_q2.clear();
      madePlots = false;  
    }





} // namespace radmat


#undef  USE_OMP_RADAMAT_DRIVER_SOLVE_LLSQ


#endif
