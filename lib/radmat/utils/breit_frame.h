#ifndef BREIT_FRAME_H_H_GUARD
#define BREIT_FRAME_H_H_GUARD

#include "itpp/itbase.h"
#include "tensor.h"

namespace radmat
{

namespace breit
{

    /**
       @brief hold a lorentz transformation and it's inverse
    */
    struct LorentzTransform
    {
        itpp::Mat<double> lambda;
        itpp::Mat<double> lambda_inv;
    };

    /**
       @brief generate the lorentz transformation to take us to the breit frame
       @details <p^prime| A | p > -> factors* <-pz | Lambda A LambdaInv | pz>
    */
    LorentzTransform genBreitLT(const itpp::Vec<double> &pmu_prime, const itpp::Vec<double> &pmu);

    //! boost two 4 vectors to the back to back frame
    LorentzTransform genBack2BackLT(const itpp::Vec<double> &pmu_prime, const itpp::Vec<double> &pmu);

    //! rapidity to boost <m,0,0,0> to <E,0,0,pz>
    double rapidity(const double E, const double p_i);

    //! generate a boost matrix in the ith cartesian direction
    itpp::Mat<double> boost_cart(const double rapidity, const int cartesian_direction);

    //! rotate the spatial components of a 4-vector r = R*o -- proper rotations
    itpp::Mat<double> rodRot4(const itpp::Vec<double> &orig, const itpp::Vec<double> &rotated);

    //! gmunu +1 time component
    itpp::Mat<double> gmunu(void);

    //! rotated = R*orig
    itpp::Mat<double> rodRotMat(const itpp::Vec<double> &orig,
                                const itpp::Vec<double> &rotated);

    //! generate the composition transformation that makes a helicity ket
    itpp::Mat<double> genH_p(const itpp::Vec<double> &pmu);

    //! generate the composition transformation that makes a helicity bra
    itpp::Mat<double> genH_p_inv(const itpp::Vec<double> &pmu);

    LorentzTransform genRotMat(const LorentzTransform &lt,
                               const itpp::Vec<double> &pmu_prim,
                               const itpp::Vec<double> &pmu);

    LorentzTransform enforceUnique(const LorentzTransform &in);

    itpp::Mat<double> R4_z(const double theta);

    // assumes a 4 vector and gets rotation for compontents orig[1-3]
    //! for the spatial compontent rotated = R*orig
    itpp::Mat<double> rodRotMat(const Tensor<double, 1> &orig,
                                const Tensor<double, 1> &rotated);

    //! return rotated = R*p_z -- default to positve z axis
    itpp::Mat<double> rodRotMat(const Tensor<double, 1> &rotated,
                                const bool plus_z = true);
} // breit

} // radmat
 
#endif
