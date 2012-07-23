#ifndef RADMAT_FIT_FORMS_H_H_GUARD
#define RADMAT_FIT_FORMS_H_H_GUARD

#include "jackFitter/jackknife_fitter.h"
#include <string>
#include <vector>

namespace radmat
{

// upon further consideration I have decided that this function 
// is absolute lunacy and I'm going home
  struct ConstPlusConstTimesTwoExp : public FitFunction
  {

    public:
      ConstPlusConstTimesTwoExp(const double tf, const double ti)
        : FitFunction(4) , m_tf(tf) , m_ti(ti) 
      {
        setParName(0,"FF");
        setParName(1,"const2");
        setParName(2,"deltamf");
        setParName(3,"deltami");
      }

      inline double operator()(const std::vector<double> &pars, double t)
      {
        return pars[0] + pars[1]*exp(pars[2]*(m_tf - t))*exp(pars[3]*(t-m_ti));
      }


      std::string getFitType(void) const {return std::string("hackey 3point fit");}

    private : 
      ConstPlusConstTimesTwoExp(void) ; // hide ctor since we need to hack jackknife to get it to work

      double m_tf;
      double m_ti;
  };

}

#endif
