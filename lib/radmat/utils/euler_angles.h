#ifndef EULER_ANGLES_H_H_GUARD
#define EULER_ANGLES_H_H_GUARD

#include "itpp/itbase.h"


namespace radmat
{
struct eAngles
{
    eAngles(void)
        : alpha(0.) , beta(0.) , gamma(0.)
    {   }

    eAngles(const eAngles &o)
        : alpha(o.alpha) , beta(o.beta) , gamma(o.gamma)
    {   }

    double alpha;
    double beta;
    double gamma;
};

// for whatever our default convention ends up as when I talk to CT
eAngles genEulerAngles(const itpp::Mat<double> &R);
eAngles genEulerAngles4(const itpp::Mat<double> &R);

// get the convention depended euler angles from a rotation matrix
eAngles genZXZ(const itpp::Mat<double> &R);  //! standard convention
eAngles genZYZ(const itpp::Mat<double> &R);  //! another convention
eAngles genXYZ(const itpp::Mat<double> &R);  //! pitch-roll-yaw convention

// do the same but pull out the lower 3x3 spatial bit of a 4x4 lorentz transformation
eAngles genZXZ4(const itpp::Mat<double> &R4);
eAngles genZYZ4(const itpp::Mat<double> &R4);
eAngles genXYZ4(const itpp::Mat<double> &R4);

}
#endif
