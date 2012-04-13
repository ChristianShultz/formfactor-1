// kinfac_key.cc -
//
// Monday, April  9 2012
//

#include "kinfac_key.h"
#include "polarisation/polarisation_factory.h"
#include "polarisation/polarisation_factory_keys.h"
#include "tensor/tensorbase.h"
#include "ensem/ensem.h"

#include <vector>
#include "omp.h"

using namespace polarisation;
using namespace kinfac;
using namespace ENSEM;

void kinfac::registerPolarisationTensors(const KinKey &key)
{

    // wrap the factory registration in a critical block
    #pragma omp critical
    {

        KinKey_p srcKey = key.src;
        KinKey_p sinkKey = key.sink;

        if(srcKey.J != 0)
            {
                int sz = srcKey.E.size();
                tensor::TensorImplBase *dum;
                tensor::Tensor<double, 1> pmu;
                idx_t J = srcKey.helicity;
                short helicity = srcKey.helicity;
                pmu.create(std::vector<idx_t>(1, 4));
                pmu[1] = srcKey.mom[0];
                pmu[2] = srcKey.mom[1];
                pmu[3] = srcKey.mom[1];
                pFac factory;
                EnsemReal E = rescaleEnsemDown(srcKey.E);

                for(int cfg = 0; cfg < sz; ++cfg)
                    {
		      pmu[0] = toDouble(E.elem(cfg));
                        dum = factory(pFacKey(pmu, J, helicity));
                    }
            }

        if(sinkKey.J != 0)
            {
                int sz = sinkKey.E.size();
                tensor::TensorImplBase *dum;
                tensor::Tensor<double, 1> pmu;
                idx_t J = sinkKey.helicity;
                short helicity = sinkKey.helicity;
                pmu.create(std::vector<idx_t>(1, 4));
                pmu[1] = sinkKey.mom[0];
                pmu[2] = sinkKey.mom[1];
                pmu[3] = sinkKey.mom[1];
                pFac factory;
                EnsemReal E = rescaleEnsemDown(sinkKey.E);

                for(int cfg = 0; cfg < sz; ++cfg)
                    {
		      pmu[0] = toDouble(E.elem(cfg));
                        dum = factory(pFacKey(pmu, J, helicity));
                    }
            }

    } // end omp critical section

}
