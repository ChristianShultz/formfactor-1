#ifndef LATTICE_MULTI_DATA_OBJECT_H
#define LATTICE_MULTI_DATA_OBJECT_H

#include "radmat/llsq/llsq_multi_data.h"
#include "lattice_multi_data_tag.h"

namespace radmat
{

  // the definition of the linear system we are creating
  typedef LLSQMultiData<LatticeMultiDataTag,std::complex<double> > LLSQLatticeMultiData; 

} // radmat

#endif /* LATTICE_MULTI_DATA_OBJECT_H */
