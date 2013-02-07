#ifndef LLSQSOLVERS_H_H_GUARD
#define LLSQSOLVERS_H_H_GUARD


#include "llsq_gen_system.h"
#include <string>
#include <complex>
#include "adat/singleton.h"
#include "adat/objfactory.h"
#include "adat/handle.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble/semble_linear_algebra.h"

namespace radmat
{

  // factory
  typedef Util::SingletonHolder<
    Util::ObjectFactory<LLSQBaseSolver_t<std::complex<double> >,
    std::string,
    void,
    LLSQBaseSolver_t<std::complex<double> >* (*)(),
    Util::StringFactoryError > >
      TheLLSQSolverFactory;

  namespace LLSQSolverFactoryEnv
  {
    bool registerAll();
    ADAT::Handle<LLSQBaseSolver_t<std::complex<double> > > callFactory(const std::string &solnID);
  }


  /////////////////////////////////////////////////////////////////////////////

  // solver functors

  // solve a square system via LU decomposition to get the inverse
  // -- this is junk and is only here as a placeholder for testing
  //    we should be using some kind of fancy SVD on non-square 
  //    overdetermined linear systems so don't use this one!!!
  template<typename T> 
    struct LLSQSolverLU_t : public LLSQBaseSolver_t<T>
  {
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h & input, const int t_ins) const
    {
      print_llsq_system(input, t_ins); 
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);
      POW2_ASSERT(input->m_KFacs.getN() == input->m_KFacs.getM()); // square linear system
      SEMBLE::SembleMatrix<T> Kinv;
      SEMBLE::inv(input->m_KFacs,Kinv);
      foo->m_FF = Kinv*(input->m_MatElems);

      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;
    }

    std::string echo(void) const {return std::string("LLSQSolverLU_t");}
  };

  // solve A*x = b for a non_square A via
  // Adag*A*x = Adag*b -- makes a square linear system
  // 
  // TO DO: write a version of this that does the whole nonsense with the residual?
  //        our system shouldn't have null space so maybe this is a waste of time..
  //

  template<typename T> 
    struct LLSQSolverSVDMakeSquare_t : public LLSQBaseSolver_t<T>
  {
    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h &input, const int t_ins) const
    {
      print_llsq_system(input, t_ins);
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);                                              // check pointer alloc 
      POW2_ASSERT(input->m_KFacs.getN() >= input->m_KFacs.getM());   // check LLSQ is not underdetermined -- it would still work
      SEMBLE::SembleMatrix<T> KinvDag_Kinv = SEMBLE::adj(input->m_KFacs)*(input->m_KFacs),U,V;
      SEMBLE::SembleVector<double> s;
      std::string svd_log = SEMBLE::svd(KinvDag_Kinv,U,s,V);
      SEMBLE::pseudoInvert(s,s.getN(),true); // s -> 1/s
      foo->m_FF = V * (SEMBLE::diagAsym<T,double>(s) )* (SEMBLE::adj(U) ) * SEMBLE::adj(input->m_KFacs) * (input->m_MatElems);

      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;

    };

    std::string echo(void) const {return std::string("LLSQSolverSVDMakeSquare_t");}

  };


  // solve A*x = b for non-square A via svd decomp of a non-square matrix
  // to do -- also return residual of the solution 

  template<typename T>
    struct LLSQSolverSVDNonSquare_t : public LLSQBaseSolver_t<T>
  {

    typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h;
    typedef typename LLSQBaseSolver_t<T>::LLSQInputType_h LLSQInputType_h;

    LLSQRetTypeBase_h operator()(const LLSQInputType_h &input, const int t_ins) const
    {
      print_llsq_system(input, t_ins);
      LLSQRetTypeBase_t<T> *foo = new LLSQRetTypeBase_t<T>();
      POW2_ASSERT(foo);                                              // check pointer alloc 
      POW2_ASSERT(input->m_KFacs.getN() >= input->m_KFacs.getM());   // check LLSQ is not underdetermined -- it would still work

      SEMBLE::SembleMatrix<T> K = input->m_KFacs;
      SEMBLE::SembleMatrix<T> U,V,Smul;
      SEMBLE::SembleMatrix<double> Sinv;
      SEMBLE::SembleVector<double> s;

      // dump the log?
      std::string svd_log = SEMBLE::svdNonSquare(K,U,s,V); 

      // do we want to cut in the future?
      SEMBLE::pseudoInvert(s,s.getN(),true); // s -> 1/s

      Sinv.reDim(s.getB(),V.getN(),U.getN());
      const int bd = (V.getN() < U.getN()) ? V.getN() : U.getN();
      Sinv.zeros();

      for(int diag = 0; diag < bd; ++diag)
        Sinv.loadEnsemElement(diag,diag,s.getEnsemElement(diag)); 

      Smul = SEMBLE::recast<T,double>(Sinv);
    
      foo->m_FF = V * Smul * SEMBLE::adj(U)  *  (input->m_MatElems);

      LLSQRetTypeBase_h ret(foo);
      print_llsq_soln(input,ret,t_ins);

      return ret;

    };

    std::string echo(void) const {return std::string("LLSQSolverSVDNonSquare_t");}


  };




}

#endif
