#ifndef AXIS_PLOTTER_H_H_GUARD
#define AXIS_PLOTTER_H_H_GUARD


#include "radmat/load_data/three_point.h"
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

  };


  template<> void MakeAxisPlots::plot(const ThreePointCorrelator<double> &C3pt,
				      const std::string &path, 
				      const std::string &fname) const;

  template<> void MakeAxisPlots::plot(const ThreePointCorrelator<std::complex<double> > &C3pt,
				      const std::string &path,
				      const std::string &fname) const;


} // namespace radmat


#endif
