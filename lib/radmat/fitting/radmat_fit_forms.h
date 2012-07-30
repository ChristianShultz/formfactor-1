#ifndef RADMAT_FIT_FORMS_H_H_GUARD
#define RADMAT_FIT_FORMS_H_H_GUARD

#include "jackFitter/jackknife_fitter.h"
#include <string>
#include <vector>

namespace radmat
{


  struct Constant: public FitFunction
  {
    Constant(void) 
      : FitFunction(1)
    {
      setParName(0,"FF");
    }

    inline double operator()(const std::vector<double> &par, double t) const
    {
      return par[0];
    }

    std::string getFitType(void) const {return std::string("Constant");}

  };


  struct ConstantTimesTwoExp : public FitFunction
  {
    ConstantTimesTwoExp (const double tf, const double ti)
      : FitFunction(3) , m_tf(tf) , m_ti(ti) 
    {
      setParName(0,"FF");
      setParName(1,"deltam_f");
      setParName(2,"deltam_i");
    }

    inline double operator()(const std::vector<double> &par, double t) const
    {
      return par[0]*exp(par[1]*(m_tf -t))*exp(par[2]*(t-m_ti));
    }

    std::string getFitType(void) const {return std::string("ConstTimesTwoExp");}

    private:
    double m_tf;
    double m_ti;
  };





} // namespace radmat

#endif
