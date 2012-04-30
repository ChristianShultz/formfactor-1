#ifndef TEST_tensor_H_H_GUARD
#define TEST_tensor_H_H_GUARD


#include "itpp/itbase.h"
#include "tensor.h"
#include "test_utils.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"

#define TEST_NARRAY_BOUND_T_GUARD 5


namespace radmat
{

    template<typename T>
    int test_tensor(void)
    {
        Tensor<T, 4> dumm4;
        Tensor<T, 4> dum4;

        std::vector<idx_t> v(4, TEST_NARRAY_BOUND_T_GUARD);

        // check creation
        POW2_ASSERT(dum4.create(&v[0]));

        // test access by setting
        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    for(idx_t l = 0; l < TEST_NARRAY_BOUND_T_GUARD; ++l)
                        dum4[i][j][k][l] = i + i * (1. - j) + k / (l + 1.) + 1. / (i + 1.);


        // test access gave what it should
        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    for(idx_t l = 0; l < TEST_NARRAY_BOUND_T_GUARD; ++l)
                        POW2_ASSERT(dum4[i][j][k][l] == i + i * (1. - j) + k / (l + 1.) + 1. / (i + 1.));

        // test assignment operator
        dumm4 = dum4;

        // test assignment operator worked properly
        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    for(idx_t l = 0; l < TEST_NARRAY_BOUND_T_GUARD; ++l)
                        POW2_ASSERT(dum4[i][j][k][l] ==  dumm4[i][j][k][l]);


        // test a 2-d contraction against matrix multiplication w/ itpp
        itpp::Mat<T> iMat1 = randu<T>(TEST_NARRAY_BOUND_T_GUARD, TEST_NARRAY_BOUND_T_GUARD), iMat2, iMat3;
        Tensor<T, 2> ctrack1, ctrack2, ctrack3;
        std::vector<idx_t> cv(2, TEST_NARRAY_BOUND_T_GUARD);
        ctrack1.create(&cv[0]);

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                ctrack1[i][j] = iMat1(i, j);

        ctrack2 = ctrack1;
        ctrack2.setIndicies(std::vector<bool>(2, false));
        ctrack3 = contract(ctrack1, ctrack2, 1, 0);
        iMat2 = iMat1;
        iMat3 = iMat1 * iMat2;

        // they apparently multiply differently, there are a few (4 for TEST_NARRAY_BOUND_T_GUARD = 5) errors of
        // order(10^-16), this will filter those
        double precision = 1e-15;

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                POW2_ASSERT(fabs(iMat3(i, j) - ctrack3[i][j]) < precision);


        // test a 4 contracting into a 4
        Tensor<T, 4> p4, pp4;
        std::vector<idx_t> v4(4, TEST_NARRAY_BOUND_T_GUARD);
        p4.create(&v4[0]);
        T p[TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD];

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    for(idx_t l = 0; l < TEST_NARRAY_BOUND_T_GUARD; ++l)
                        {
                            p[i][j][k][l] = i + i * (1. - j) + k / (l + 1.) + 1. / (i + 1.);
                            p4[i][j][k][l] = p[i][j][k][l];
                        }

        T p2[TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD]
        [TEST_NARRAY_BOUND_T_GUARD];

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    for(idx_t l = 0; l < TEST_NARRAY_BOUND_T_GUARD; ++l)
                        for(idx_t m = 0; m < TEST_NARRAY_BOUND_T_GUARD; ++m)
                            for(idx_t n = 0; n < TEST_NARRAY_BOUND_T_GUARD; ++n)
                                {
                                    p2[i][j][k][l][m][n] = 0.;

                                    for(idx_t c = 0; c < TEST_NARRAY_BOUND_T_GUARD; ++c)
                                        p2[i][j][k][l][m][n] += p[i][j][c][k] * p[l][c][m][n];
                                }

        pp4 = p4;
        pp4.flipAllIndicies();
        Tensor<T, 6> p6 = contract(pp4, p4, 2, 1);

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    for(idx_t l = 0; l < TEST_NARRAY_BOUND_T_GUARD; ++l)
                        for(idx_t m = 0; m < TEST_NARRAY_BOUND_T_GUARD; ++m)
                            for(idx_t n = 0; n < TEST_NARRAY_BOUND_T_GUARD; ++n)
                                POW2_ASSERT(fabs(p2[i][j][k][l][m][n] - p6[i][j][k][l][m][n]) < precision);


        // test the algebraic operations
        Tensor<T, 3> o3base, omul, odiv, oadd, osub, ominus;
        std::vector<idx_t> o3dim(3, TEST_NARRAY_BOUND_T_GUARD);
        o3base.create(&o3dim[0]);

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    o3base[i][j][k] = i + i * (1. - j) + k + 1. / (1. + i);

        omul = o3base * 2.;
        odiv = o3base / 2.;
        oadd = o3base + o3base;
        osub = o3base - o3base;
        ominus = -o3base;

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    {
                        POW2_ASSERT(omul[i][j][k] == 2.* o3base[i][j][k]);
                        POW2_ASSERT(odiv[i][j][k] == .5 * o3base[i][j][k]);
                        POW2_ASSERT(oadd[i][j][k] == o3base[i][j][k] + o3base[i][j][k]);
                        POW2_ASSERT(osub[i][j][k] == 0.);
                        POW2_ASSERT(ominus[i][j][k] == -o3base[i][j][k]);
                    }

        // test the array slicing
        Tensor<T, 2> slice;

        for(int cut = 0; cut < TEST_NARRAY_BOUND_T_GUARD; ++cut)
            {
                slice = o3base.slice(1, cut);

                for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
                    for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                        POW2_ASSERT(slice[i][j] == o3base[i][cut][j]);
            }

        // test the array contraction for metrics
        Tensor<T, 2> metric;
        std::vector<idx_t> mdim(2, TEST_NARRAY_BOUND_T_GUARD);
        metric.create(&mdim[0]);
        metric.flipAllIndicies();

        // create a metric of -I_(TEST_NARRAY_BOUND_T_GUARD-1) and then check
        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                if(i == j)
                    metric[i][j] = -1.;
                else
                    metric[i][j] = 0.;

        metric[TEST_NARRAY_BOUND_T_GUARD - 1][TEST_NARRAY_BOUND_T_GUARD - 1] = 0;

        Tensor<T, 3> bar;
        bar = applyMetric(o3base, metric, 1);

        for(idx_t i = 0; i < TEST_NARRAY_BOUND_T_GUARD; ++i)
            for(idx_t j = 0; j < TEST_NARRAY_BOUND_T_GUARD; ++j)
                for(idx_t k = 0; k < TEST_NARRAY_BOUND_T_GUARD; ++k)
                    if(j != TEST_NARRAY_BOUND_T_GUARD - 1)
                        POW2_ASSERT(bar[i][j][k] == -o3base[i][j][k]);
                    else
                        POW2_ASSERT(bar[i][j][k] == 0.);

        std::cout << "All " << Stringify<T>() << " tensor tests successfull" << std::endl;

        return 0;
    }

}

#undef TEST_NARRAY_BOUND_T_GUARD


#endif
