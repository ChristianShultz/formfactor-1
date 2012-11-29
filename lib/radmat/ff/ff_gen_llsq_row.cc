// ff_gen_llsq_row.cc -
//
// Monday, July  9 2012
//

#include"ff_gen_llsq_row.h"
#include "radmat/utils/tensor.h"
#include "semble/semble_vector.h"
#include "radmat/utils/pow2assert.h"

using namespace SEMBLE;

namespace radmat
{

  Tensor<double,1> pPlus(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i)
  {
    return p_f + p_i;
  }

  Tensor<double,1> pMinus(const Tensor<double,1> &p_f, const Tensor<double,1> &p_i)
  {
    return p_f - p_i;
  }

 SemblePInv makeMomInvariants(const EnsemReal &E_f, const EnsemReal &E_i, const Array<int> &p_f, const Array<int> &p_i, const double factor)
 {
   SembleVector<double> final(E_f.size(),4);
   SembleVector<double> initial(E_i.size(),4);

   POW2_ASSERT(final.getB() == initial.getB());

   int ncfg = final.getB();

   final.loadEnsemElement(0,E_f);
   initial.loadEnsemElement(0,E_i);


   for(int bin = 0; bin < ncfg; bin ++)
     for(int lorentz = 0; lorentz < 3; lorentz ++)
       {
	 final.setElement(bin,lorentz+1,factor*double(p_f[lorentz]));
	 initial.setElement(bin,lorentz+1,factor*double(p_i[lorentz]));
       }

    SemblePInv p;
    p.pi() = initial;
    p.pf() = final; 

    return p;
 };

} // close radmat
