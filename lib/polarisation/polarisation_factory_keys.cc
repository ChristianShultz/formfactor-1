// polarisation_factory_keys.cc -
//
// Wednesday, March 28 2012
//

#include "polarisation_factory_keys.h"
#include <iostream>

using namespace polarisation;
using namespace tensor;

pFacKey::pFacKey(Tensor<double,1> &_pmu, idx_t _J, short _helicity)
  : pmu(_pmu), J(_J) , helicity(_helicity)
{    }


pFacKey& pFacKey::operator=(const pFacKey &o)
{
  if(this != &o)
    {
      pmu = o.pmu;
      J = o.J;
      helicity = o.helicity;
    }
  return *this;
}

double& pFacKey::operator[](const idx_t i)
{
  return pmu[i];
}

const double& pFacKey::operator[](const idx_t i) const
{
  return pmu[i];
}

// stupid weak ordering
bool pFacKeyComp::operator()(const pFacKey &lhs, const pFacKey &rhs)
{
  return lhs < rhs;
}

bool polarisation::operator<(const pFacKey &lhs, const pFacKey &rhs)
{
  if(lhs.pmu != rhs.pmu)
    {
      idx_t vec_len = lhs.pmu.getDim(0);

      for(idx_t i = 0; i < vec_len; ++i) 
	if(lhs[i] != rhs[i])
	  return lhs[i] < rhs[i];
    }

  if(lhs.J != rhs.J)
    return lhs.J < rhs.J;

  return lhs.helicity < rhs.helicity;
}

bool polarisation::operator>(const pFacKey &lhs, const pFacKey &rhs)
{
  if(lhs.pmu != rhs.pmu)
    {
      idx_t vec_len = lhs.pmu.getDim(0);

      for(idx_t i = 0; i < vec_len; ++i) 
	if(lhs.pmu[i] != rhs.pmu[i])
	  return lhs.pmu[i] > rhs.pmu[i];
    }

  if(lhs.J != rhs.J)
    return lhs.J > rhs.J;

  return lhs.helicity > rhs.helicity;
}

bool polarisation::operator==(const pFacKey &lhs, const pFacKey &rhs)
{
  return !(lhs < rhs) && !(lhs > rhs);
}

bool polarisation::operator!=(const pFacKey &lhs, const pFacKey &rhs)
{
  return !(lhs == rhs);
}


std::ostream& polarisation::operator<<(std::ostream &o, const pFacKey &k )
{
  o << "J = " << k.J << "\n";
  o << "helicity = " << k.helicity << "\n";
  o << "pmu = " << k.pmu;
  return o;
}
