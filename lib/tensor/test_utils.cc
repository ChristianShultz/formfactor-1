// test_utils.cc -
//
// Monday, March 26 2012
//

#include "test_utils.h"
#include "itpp/itbase.h"
#include <complex>


// specializations and an overload for tensor::fabs since on my computer it didnt pick it up from
// the std namespace, why?
namespace tensor
{

    template<>
    itpp::Mat<double> randu(const int a, const int b)
    {
        return itpp::randu(a, b);
    }

    template<>
    itpp::Mat<std::complex<double> > randu(const int a, const int b)
    {
        return itpp::randn_c(a, b);
    }

    double fabs(std::complex<double> c)
    {
        return abs(c);
    }

}
