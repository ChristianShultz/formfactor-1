// euler_angles.cc -
//
// Tuesday, April 24 2012
//

#include "euler_angles.h"
#include "pow2assert.h"
#include <math.h>

namespace radmat
{

  eAngles genEulerAngles(const itpp::Mat<double> &R)
  {
    return genZXZ(R);
  }

  eAngles genEulerAngles4(const itpp::Mat<double> &R)
  {
    return genZXZ4(R);
  }

  eAngles genZXZ(const itpp::Mat<double> &R)
  {
    POW2_ASSERT((R.rows() == R.cols()) && (R.rows() == 3));

    eAngles ret;


    ret.alpha = -atan(R(2, 0) / R(2, 1));
    ret.beta = acos(R(2, 2));
    ret.gamma = atan(R(0, 2) / R(1, 2));

    return ret;
  }

  eAngles genZYZ(const itpp::Mat<double> &R)
  {
    eAngles ret = genZXZ(R);
    const static double pi = acos(-1.);
    ret.alpha -= pi / 2.;
    ret.gamma += pi / 2.;

    return ret;
  }


  eAngles genXYZ(const itpp::Mat<double> &R)
  {
    eAngles ret;

    ret.alpha = atan(R(0, 1) / R(0, 0));
    ret.beta = asin(-R(0, 3));
    ret.gamma = atan(R(1, 2) / R(2, 2));

    return ret;
  }


  eAngles genZXZ4(const itpp::Mat<double> &foo)
  {
    itpp::Mat<double> foobar(foo);
    foobar.del_row(0);
    foobar.del_col(0);

    return genZXZ(foobar);
  }

  eAngles genZYZ4(const itpp::Mat<double> &foo)
  {
    itpp::Mat<double> foobar(foo);
    foobar.del_row(0);
    foobar.del_col(0);

    return genZYZ(foobar);
  }

  eAngles genXYZ4(const itpp::Mat<double> &foo)
  {
    itpp::Mat<double> foobar(foo);
    foobar.del_row(0);
    foobar.del_col(0);

    return genXYZ(foobar);
  }

}// radmat
