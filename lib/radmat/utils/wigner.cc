// wigner.cc -
//
// Tuesday, April 24 2012
//

#include"wigner.h"
#include "breit_frame.h"
#include "hadron/clebsch.h"
#include "semble/semble_meta.h"
#include "euler_angles.h"
#include "tensor.h"
#include <math.h>
#include <complex>
#include "pow2assert.h"


using namespace std;

namespace radmat
{

namespace wigner
{

    namespace
    {
        // adat/lib/hadron/irrep_utils.cc -- looks like CT chose D( R(alpha,beta,0))
        // in his definition of continuum helicity operators
        // changing this will apply whatever convention we want
        void wignerApplyConvention_(eAngles &a)
        {
            a.gamma = 0.;
        }
    } // end anonymous namespace

    std::complex<double> Wigner_D(const int J, const int M, const int N, const eAngles &eul_angles)
    {
        // always apply the rotation conventions that christopher used
        eAngles my_eul(eul_angles);
        wignerApplyConvention_(my_eul);
        return SEMBLE::toScalar(Hadron::Wigner_D(2 * J, 2 * M, 2 * N, my_eul.alpha, my_eul.beta, my_eul.gamma));
    }

    std::complex<double> Wigner_D_noConvention(const int J, 
					       const int M, 
					       const int N, 
					       const eAngles &eul)
    {
        return SEMBLE::toScalar(Hadron::Wigner_D(2 * J, 2 * M, 2 * N, eul.alpha, eul.beta, eul.gamma));
    }

    std::complex<double> Wigner_D(const int J, 
				  const int M, 
				  const int N, 
				  const Tensor<double, 1> &pmu, bool pz)
    {
        return Wigner_D(J, M, N, breit::rodRotMat(pmu, pz));
    }

    std::complex<double> Wigner_D(const int J, 
				  const int M, 
				  const int N, 
				  const itpp::Mat<double> &R)
    {
        if(R.rows() == 4)
            return Wigner_D(J, M, N, genEulerAngles4(R));

        if(R.rows() == 3)
            return Wigner_D(J, M, N, genEulerAngles(R));

        // shouldn't get here unless something goes wrong, then we can backtrace the exit in a debugger to
        // find the offending section of code
        POW2_ASSERT(false); // force a bad instruction causing an abort

	return std::complex<double>(10^100,10^100); // force the compiler to shut up about the warning
    }


  double Wigner_3j(const int j1, 
		   const int j2, 
		   const int J, 
		   const int m1, 
		   const int m2, 
		   const int m)
    {
        int phase = (j1 - j2 - m) % 2;
        double factor = 1. / sqrt(2.*J + 1) * ((-1) ^ phase);
	//clebsch takes 2 times spin
        return (factor * Hadron::clebsch(2 * j1, 2 * m1, 2 * j2, 2 * m2, 2 * J, -2 * m) ); 
    }

}// close wigner namespace


} // radmat
