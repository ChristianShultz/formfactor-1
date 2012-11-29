#ifndef AXIS_PLOTTER_H_H_GUARD
#define AXIS_PLOTTER_H_H_GUARD


#include "radmat/load_data/three_point.h"
#include "radmat/fake_data/fake_3pt_function.h"
#include <string>
#include <complex>



namespace radmat
{

  struct
    MakeAxisPlots
    {
      template<typename T>
        void plot(const ThreePointCorrelator<T> &C3pt,
            const std::string &path,
            const std::string &fname) const;


      template<typename T>
        void plot(const Fake3ptCorr<T> &) const; 

    };


  template<> void MakeAxisPlots::plot(const ThreePointCorrelator<double> &C3pt,
      const std::string &path, 
      const std::string &fname) const;

  template<> void MakeAxisPlots::plot(const ThreePointCorrelator<std::complex<double> > &C3pt,
      const std::string &path,
      const std::string &fname) const;

  template<> void MakeAxisPlots::plot(const Fake3ptCorr<double> &)const;

  template<> void MakeAxisPlots::plot(const Fake3ptCorr<std::complex<double> >&) const;





} // namespace radmat


#endif
