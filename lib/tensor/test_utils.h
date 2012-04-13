#ifndef TEST_UTILS_H_H_GUARD
#define TEST_UTILS_H_H_GUARD

#include "itpp/itbase.h"

#include <complex>

namespace tensor
{
    template<typename T>
    itpp::Mat<T> randu(const int, const int);

    double fabs(std::complex<double> c);
}
#endif
