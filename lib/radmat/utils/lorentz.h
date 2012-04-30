#ifndef LORENTZ_H_H_GUARD
#define LORENTZ_H_H_GUARD

#include "itpp/itbase.h"
#include "semble/semble_matrix.h"
#include "euler_angles.h"
#include <utility>

/**
   @file lorentz.h
   @brief some lorentz transformation linear algebra
 */

// NB all of this assumes we are getting 4 vectors with raised indicies
namespace radmat
{

namespace lorentz
{

    /**
       @brief hold a transformation and it's inverse
     */
    struct LorentzTransform
    {
        itpp::Mat<double> Lambda;
        itpp::Mat<double> LambdaInv;
    };


    /**
       @brief hold an ensemble of lorentz transformations and do jackknife
     */
    struct EnsembleLT
    {
        EnsembleLT(const int nbins)
        {
            Lambda.reDim(nbins, 4, 4);
            LambdaInv = Lambda;
        }

        void reDim(const int nbins)
        {
            Lambda.reDim(nbins, 4, 4);
            LambdaInv = Lambda;
        }

        void scaleUp(void)
        {
            Lambda.rescaleEnsemUp();
            LambdaInv.rescaleEnsemUp();
        }

        void scaleDown(void)
        {
            Lambda.rescaleEnsemDown();
            LambdaInv.rescaleEnsemDown();
        }

        void add(const lorentz::LorentzTransform &LT, const int bin)
        {
            Lambda[bin] = LT.Lambda;
            LambdaInv[bin] = LT.LambdaInv;
        }

        lorentz::LorentzTransform get(const int bin)
        {
            lorentz::LorentzTransform foo;
            foo.Lambda = Lambda[bin];
            foo.LambdaInv = LambdaInv[bin];
            return foo;
        }

        SEMBLE::SembleMatrix<double> Lambda;
        SEMBLE::SembleMatrix<double> LambdaInv;
    };

    /**
       @brief generate the lorentz transformation to take us to the breit frame
       @details <p^prime| A | p > -> factors* <-pz | Lambda A LambdaInv | pz>
     */
    LorentzTransform genBreitLT(const itpp::Vec<double> &p_prime4, const itpp::Vec<double> &p4);

    //! rapidity to boost <m,0,0,0> to <E,0,0,pz>
    double rapidity(const double E, const double pz);

    //! return a rotation matrix to rotate the spatial components of p4 to z-axis (negative z-axis if false)
    itpp::Mat<double> rotateToZAxis(const itpp::Vec<double> &p4, const bool plus = true);

    //! return the boost matrix that takes <E,0,0,pz> to <m,0,0,0>
    itpp::Mat<double> boostZ(const double rapidity);

    //! boost in the ith direction
    itpp::Mat<double> boostI(const double rapidity, int i);

    //! x-convention is used
    itpp::Mat<double> rotMatrix4(const itpp::Vec<double> &euler_angles);

    //! generate the Euler Angles to rotate spatial components of v4 to the z-axis
    //! using the x-convention
    itpp::Vec<double> genEulerAnglesToZAxis(const itpp::Vec<double> &v4);

    //! generate the Euler Angles to rotate spatial components of v4 to the negative z-axis
    //! using the x-convention
    itpp::Vec<double> genEulerAnglesToNegativeZAxis(const itpp::Vec<double> &v4);

    //! gmunu +1 temporal
    itpp::Mat<double> gmunu(void);

    //! return the pair of euler angles according to Durand
    std::pair<eAngles, eAngles> eulerAngles(const LorentzTransform &lt,
                                            const itpp::Vec<double> &pp4,
                                            const itpp::Vec<double> &p4);

} //lorentz

} //radmat
#endif
