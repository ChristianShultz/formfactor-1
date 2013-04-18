/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name :

* Purpose :

* Creation Date : 14-02-2013

* Last Modified : Thu Feb 14 09:29:16 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/


#include "fake_3pt_function_aux.h"
#include <math.h>



namespace radmat
{

  double mom_factor(const double xi, const int L_s) 
  {
    return 2.*acos(-1.)/xi/double(L_s);
  }

}
