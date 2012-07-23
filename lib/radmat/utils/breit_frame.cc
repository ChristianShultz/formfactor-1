// breit_frame.cc -
//
// Wednesday, April 25 2012
//

#include "breit_frame.h"
#include "itpp/itbase.h"
#include "tensor.h"

using namespace std;
using namespace itpp;

namespace radmat
{

namespace breit
{

    namespace
    {
        double precision = 1e-14;

        double breit_get_mass_(const Vec<double> &pmu)
        {
            return sqrt(pmu * (gmunu() * pmu));
        }

        double breit_get_modp_(const Vec<double> &pmu)
        {
            return sqrt(pmu[1] * pmu[1] + pmu[2] * pmu[2] + pmu[3] * pmu[3]);
        }

        bool breit_within_tolerance_(const itpp::Vec<double> &a, const itpp::Vec<double> &b, const double tolerance)
        {
            if(a.size() != b.size())
                return false;

            for(int i = 0; i < a.size(); ++i)
                if(fabs(a[i] - b[i]) > tolerance)
                    return false;

            return true;
        }

        // the rotation may not be improper if there is a zero element, deal with this case
        itpp::Mat<double> breit_rod_rot_opposite_(const itpp::Vec<double> &o,
                const itpp::Vec<double> &r,
                const double tol)
        {
            // find a zero element if it exists so we can rotate about that axis by pi
            int zidx(10);

            for(int i = 0; i < 3; ++i)
                if(fabs(o[i]) < tol)
                    zidx = i;

            itpp::Vec<double> I3(3);
            I3.ones();
            itpp::Mat<double> mI(-diag(I3));

            // didn't find one, return an improper rotation
            if(zidx == 10)
                {
                    cout << __func__ << " warning: improper rotation performed" << endl;
                    return mI;
                }

            mI(zidx, zidx) = 1;
            return mI;
        }

    } // end anonymous namespace


    double rapidity(const double E, const double pz)
    {
        return 0.5 * log((E - pz) / (E + pz));
    }


    itpp::Mat<double> boost_cart(const double xi, const int i)
    {
        POW2_ASSERT((i < 4) && (i > 0));
        itpp::Mat<double> K(4, 4);
        K.zeros();
        itpp::Vec<double> I4(4);
        I4.ones();

        K(i, 0) = 1;
        K(0, i) = 1;

        return (
                   itpp::diag(I4)
                   - K * K
                   - K * sinh(xi)
                   + K * K * cosh(xi)
               );
    }

    LorentzTransform genBack2BackLT(const itpp::Vec<double> &pmu_prime, const itpp::Vec<double> &pmu)
    {
        POW2_ASSERT((pmu.size() == pmu_prime.size()) && (pmu.size() == 4));

        itpp::Vec<double> pplus = pmu + pmu_prime;

        itpp::Mat<double> B1, B2, B3, Bi1, Bi2, Bi3;
        double r;

        r = rapidity(pplus[0], pplus[1]);
        B1 = boost_cart(-r , 1);
        Bi1 = boost_cart(r , 1);
        pplus = B1 * pplus;

        r = rapidity(pplus[0], pplus[2]);
        B2 = boost_cart(-r , 2);
        Bi2 = boost_cart(r , 2);
        pplus = B2 * pplus;

        r = rapidity(pplus[0], pplus[3]);
        B3 = boost_cart(-r , 3);
        Bi3 = boost_cart(r , 3);
        pplus = B3 * pplus;

        LorentzTransform LT;

        LT.lambda = B3 * B2 * B1;
        LT.lambda_inv = Bi1 * Bi2 * Bi3;

        /*
        // DEBUG
        cout << "pplus        " << pplus << endl;
        cout << "B*pmu        " << LT.lambda * pmu << endl;
        cout << "B*pmu_prime  " << LT.lambda * pmu_prime << endl;

        // DEBUG
        cout << breit_get_mass_(LT.lambda * pmu) << endl;
        cout << breit_get_mass_(LT.lambda * pmu_prime) << endl;
        */

        return LT;
    }

    itpp::Mat<double> rodRot4(const itpp::Vec<double> &orig, const itpp::Vec<double> &rotated)
    {
        POW2_ASSERT((orig.size() == rotated.size()) && (orig.size() == 4));

        itpp::Vec<double> o(3) , r(3) ;

        for(short i = 0; i < 3; ++i)
            {
                o[i] = orig[i + 1];
                r[i] = rotated[i + 1];
            }

        itpp::Mat<double> R(4, 4) , rmat = rodRotMat(o, r);
        R.zeros();
        R(0, 0) = 1.;

        for(short i = 0; i < 3; ++i)
            for(short j = 0; j < 3; ++j)
                R(i + 1, j + 1) = rmat(i, j);

        return R;
    }

    // couldn't find it in a real book but wolfram came to the rescue
    // http://mathworld.wolfram.com/RodriguesRotationFormula.html

    // solves rotated = R*orig for the rotation matrix R where R is a proper rotation
    itpp::Mat<double> rodRotMat(const itpp::Vec<double> &orig,
                                const itpp::Vec<double> &rotated)
    {
        POW2_ASSERT((orig.size() == rotated.size()) && (orig.size() == 3));

        itpp::Vec<double> no = orig / sqrt(orig * orig);
        itpp::Vec<double> nr = rotated / sqrt(rotated * rotated);
        itpp::Vec<double> cross = itpp::cross(no, nr), I3(3);
        itpp::Mat<double>  cross3(3, 3);
        I3.ones();

        // the general formula will fail for improper rotations, ie parity, we can explicitly make
        // this work by returning -I3 in this case

        if(breit_within_tolerance_(no, -nr, precision)) // opposite vectors
            return breit_rod_rot_opposite_(no, nr, precision);

        double theta = acos(no * nr);

        // need to renormalize this guy to a unit vector
        if(cross * cross > precision)     // numerically compatible with a zero vector -- ~colinear vectors
            cross /= sqrt(cross * cross); // this should never happen with the lattice p vectors

        /*
          std::cout << "no * no = " << no * no << std::endl;
          std::cout << "nr * nr = " << nr * nr << std::endl;
          std::cout << "no * cross = " << no*cross << std::endl;
          std::cout << "nr * cross = " << nr*cross << std::endl;
          std::cout << "cross * cross = " << cross*cross << std::endl;
        */

        cross3.zeros();

        cross3(0, 1) = -cross(2);
        cross3(1, 0) = cross(2);
        cross3(0, 2) = cross(1);
        cross3(2, 0) = -cross(1);
        cross3(1, 2) = -cross(0);
        cross3(2, 1) = cross(0);

        return itpp::diag(I3) + cross3 * sin(theta) + cross3 * cross3 * (1. - cos(theta));
    }

    itpp::Mat<double> gmunu(void)
    {
        itpp::Vec<double> g(4) ;
        g.ones();
        g = - g;
        g(0) = 1;
        return diag(g);
    }


    LorentzTransform genBreitLT(const itpp::Vec<double> &pmu_prime , const itpp::Vec<double> &pmu)
    {
        LorentzTransform lt = genBack2BackLT(pmu_prime, pmu);
        itpp::Vec<double> pmu_b2b, pmup_b2b, z(4);

        pmu_b2b = round_to_zero(lt.lambda * pmu, precision);

        z.zeros();
        z[3] = 1.;

        itpp::Mat<double> R = rodRot4(pmu_b2b, z);

        lt.lambda = R * lt.lambda;
        lt.lambda_inv = lt.lambda_inv * transpose(R);

        // DEBUG
        // cout << __func__ << " " << lt.lambda * pmu << endl;
        // cout << __func__ << " " << lt.lambda * pmu_prime << endl;
        // cout << __func__ << "\n" << lt.lambda << endl;
        // cout << __func__ << "\n" << lt.lambda << endl;

        return enforceUnique(lt);
    }

    itpp::Mat<double> genH_p(const itpp::Vec<double> &pmu)
    {
        POW2_ASSERT(pmu.size() == 4);
        double m = breit_get_mass_(pmu);
        double p = breit_get_modp_(pmu);
        Vec<double> my_pmu(4), rest;
        my_pmu.zeros();
        itpp::Mat<double> boost = boost_cart(rapidity(pmu[0], p) , 3);
        my_pmu[0] = m;
        rest = my_pmu;
        my_pmu = boost * my_pmu;

        /*  // DEBUG
        std::cout << __func__ << " pmu " << pmu << std::endl;
        std::cout << __func__ << " B(p)*rest " << boost * rest<< std::endl;
        std::cout << __func__ << " R(p)*B(p)*rest " << rodRot4(my_pmu,pmu) *boost * rest << std::endl;
        std::cout << __func__ << " H(p) \n" << rodRot4(my_pmu,pmu)*boost << std::endl;
        std::cout << __func__ << "H_inv(p) * H(p) * rest "
        << genH_p_inv(pmu) * rodRot4(my_pmu,pmu)*boost * rest << std::endl;
        */

        return (rodRot4(my_pmu, pmu) * boost);
    }


    itpp::Mat<double> genH_p_inv(const itpp::Vec<double> &pmu)
    {
        POW2_ASSERT(pmu.size() == 4);
        double m = breit_get_mass_(pmu);
        double p = breit_get_modp_(pmu);
        Vec<double> my_pmu(4);
        my_pmu.zeros();
        itpp::Mat<double> boost = boost_cart(-rapidity(pmu[0], p), 3);
        my_pmu[0] = m;
        my_pmu = boost_cart(rapidity(pmu[0], p), 3) * my_pmu;

        return boost * transpose(rodRot4(my_pmu, pmu));
    }


    LorentzTransform genRotMat(const LorentzTransform &lt,
                               const itpp::Vec<double> &pmup,
                               const itpp::Vec<double> &pmu)
    {
        LorentzTransform LT;

        itpp::Vec<double> pbreit, ppbreit;

        pbreit = lt.lambda * pmu;
        ppbreit = lt.lambda * pmup;

        LT.lambda = genH_p_inv(pbreit) * lt.lambda * genH_p(pmu);
        LT.lambda_inv = genH_p_inv(pmup) * lt.lambda_inv * genH_p(ppbreit);

        // DEBUG
        // std::cout << __func__ << " pbreit " << pbreit << std::endl;
        // std::cout << __func__ << " ppbreit " << ppbreit << std::endl;
        // std::cout << LT.lambda << std::endl;
        // std::cout << LT.lambda_inv << std::endl;

        return LT;
    }

    LorentzTransform enforceUnique(const LorentzTransform &in)
    {
        LorentzTransform out;
        double theta = atan(in.lambda(1, 1) / in.lambda(2, 1));
        itpp::Mat<double> R = R4_z(theta);

        out.lambda = R * (in.lambda);
        out.lambda_inv = (in.lambda_inv) * transpose(R);

        return out;
    }


    itpp::Mat<double> R4_z(const double theta)
    {
        itpp::Mat<double> foo(4, 4);
        foo.zeros();
        foo(0, 0) = foo(3, 3) = 1.;
        foo(1, 1) = foo(2, 2) = cos(theta);
        foo(1, 2) = -sin(theta);
        foo(2, 1) = -foo(1, 2);

        return foo;
    }

    itpp::Mat<double> rodRotMat(const Tensor<double, 1> &orig, const Tensor<double, 1> &rot)
    {
        itpp::Vec<double> o(3), r(3);

        for(idx_t i = 0; i < 3; ++i)
            {
                o[i] = orig[i + 1];
                r[i] = rot[i + 1];
            }

        return rodRotMat(o, r);
    }

    itpp::Mat<double> rodRotMat(const Tensor<double, 1> &r, bool pz)
    {
      Tensor<double, 1> z((TensorShape<1>())[4],0.);
      z[3] = pz ? 1. : -1.;
      return rodRotMat(z, r);
    }

} // close breit namespace

} // radmat
