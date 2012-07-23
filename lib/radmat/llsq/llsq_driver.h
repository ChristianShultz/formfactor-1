#ifndef LLSQ_DRIVER_H_H_GUARD
#define LLSQ_DRIVER_H_H_GUARD


#include "llsq_q2_pack.h"
#include "llsq_solver.h"
#include "adat/handle.h"
#include <string>

namespace radmat
{


  template<typename T>
    struct LLSQDriver_t
    {
    // save some typing
    typedef typename LLSQSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQSolver_t<T>::LLSQInputType_h LLSQInputType_h;


      LLSQDriver_t(const std::string &solverID)
        : m_solverID(solverID)
      {  }


      LLSQRet_ff_Q2Pack<T> operator()(const ADAT::Handle<LLSQDataPointQ2Pack> &in) const
      {
        LLSQDataPointQ2Pack::const_iterator it;
        LLSQRet_t_Q2Pack<T> tmp;

        for(it = in->begin(); it != in->end(); it++)
          tmp.insert(it->first, m_solver(it->second, m_solverID));


        tmp.setQ2(in->Q2());

        return transformLLSQRetPack(tmp);
      }


      private:
      LLSQDriver_t(void); 
      LLSQSolver_t<T> m_solver;
      std::string m_solverID;
    };



} // radmat






#endif
