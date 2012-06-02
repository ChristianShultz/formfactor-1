#ifndef LLSQ_SOLVER_H_H_GUARD
#define LLSQ_SOLVER_H_H_GUARD


#include "LLSQ_base.h"
#include "LLSQSolvers.h"
#include "adat/handle.h"
#include <string>

namespace radmat
{

  // i hate xml ini files but we obviously need one here..
  struct LLSQIni_t
  {
    std::string solnID;
    std::string matElemID;
    double Q2PrecisionCheck; // set to -1 to avoid a check
  };


  // this basically sets up some matrix element/soltution type dependent linear system and then 
  // returns the solution
  template<typename T>
  struct LLSQSolver_t
  {
    // save some typing
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;
    typedef typename ADAT::Handle<LLSQBaseSolver_t<T> > LLSQBaseSolver_h; 
    
    // do nothing
    LLSQSolver_t(void) {};
    // another functor 
    LLSQRetTypeBase_h operator()(const typename LLSQInputType_t<T>::SemblePInvList_t & momenta,
				 const typename LLSQInputType_t<T>::LatticeMatrixElements &matElems,
				 const LLSQIni_t &iniKeys)
    {
      // if we want to mess around with what goes in we can throw a key into the LLSQIni_t and do it here
      // obviously we are going to need to dynamically cast out of the base handle if we want to send in 
      // some fancy type of inputs that are solution method dependet but this should be simple/safe to implement

      // strange that "solving" the problem is actually done in three lines
      LLSQInputType_h inputs(new LLSQInputType_t<T>(momenta,matElems,inikeys.matElemID,inikeys.Q2PrecisionCheck));
      LLSQBaseSolver_h solver = LLSQSolverFactoryEnv::callFactory(inikeys.solnID);

      return (*solver)(inputs);
    }				     
  };



}

#endif
