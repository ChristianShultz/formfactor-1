#ifndef WIGNER_H_H_GUARD
#define WIGNER_H_H_GUARD


#include <complex>
#include "semble/semble_meta.h"
#include "euler_angles.h"
#include "radmat/tensor/tensor.h"
#include "itpp/itbase.h"

namespace radmat
{
namespace wigner
{

  /** use the hadron/clebsh but allow for the different rotation conventions with an easy switch here

      note: the args are not 2*spin and this will only work for bosons

      this returns the D matrix for the rotation from pz to pmu, or -pz to pmu
  */

  // NB: in this namespace we are dealing with mesons and thus all arguments are spin, the factor of 2
  //     to line up with adat is included by the functions

  std::complex<double> Wigner_D(const int J,
				const int M,
				const int N,
				const Tensor<double, 1> &pmu,
				bool plus_z = true);

  std::complex<double> Wigner_D(const int J,
				const int M,
				const int N,
				const itpp::Mat<double> &R);

  /**
     piggyback off the routine in hadron and stick in the factors of 2, apply gamma = 0 convention
  */
  std::complex<double> Wigner_D(const int J,
				const int M,
				const int N,
				const eAngles &eulerAngles);

  //!  piggyback off the routine in hadron and stick in the factors of 2, don't apply the gamma = 0 convention
  std::complex<double> Wigner_D_noConvention(const int J,
					     const int M,
					     const int N,
					     const eAngles &eulerAngles);

  //! return the value of the wigner 3-J symbol
  double Wigner_3j(const int j1, const int j2, const int J, const int m1, const int m2, const int m);



} //wigner

} //radmat
#endif
