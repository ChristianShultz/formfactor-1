#ifndef LLSQ_GEN_SYSTEM_H_H_GUARD
#define LLSQ_GEN_SYSTEM_H_H_GUARD

#include "ff_gen_llsq_row.h"
#include "formfactor_factory.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/splash.h"
#include "semble/semble_vector.h"
#include <vector>

using namespace ADATXML;
using namespace ADATIO;

namespace radmat
{

  struct LLSQDataPoint
  {
    std::string matElemID;                      // int the J_f J_i language
    std::pair<bool,ENSEM::EnsemComplex> zero;   // lorentz index of measurements
    std::pair<bool,ENSEM::EnsemComplex> one;    // bool is if we want to use it
    std::pair<bool,ENSEM::EnsemComplex> two;
    std::pair<bool,ENSEM::EnsemComplex> three;
    Array<int> p_f;                              // the momenta used
    Array<int> p_i;
    EnsemReal E_f;                              // these are in flight energies 
    EnsemReal E_i;                              //  ie: sqrt(m*m + p*p)
  };

  /////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // deriving from this polymorphicall will inevitably 
  // give us some freedom later, return handles from the solver routines
  template<typename T>
  struct LLSQRetTypeBase_t
  {
    typedef typename SEMBLE::SembleVector<T> FFType; 
        
    LLSQRetTypeBase_t(const FFType &FFSoln)
    : m_FF(FFSoln)
    {  }

    LLSQRetTypeBase_t(void) 
    : m_FF(SEMBLE::SembleVector<T>())
    {  }

    virtual ~LLSQRetTypeBase_t(void) {}

  public: 
    FFType m_FF;
  };


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  template<typename T>
  struct LLSQInputType_t
  {
    typedef typename SEMBLE::SembleVector<T> LatticeMatrixElements;
    typedef typename SEMBLE::SembleMatrix<T> KinematicFactors;

    LLSQInputType_t(const KinematicFactors &KFacs, const LatticeMatrixElements &MatElems)
      : m_KFacs(KFacs) , m_MatElems(MatElems)
    {  }

  private:
    LLSQInputType_t(void);

  public:
    KinematicFactors m_KFacs;
    LatticeMatrixElements m_MatElems;
  };

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  template<typename T>
  LLSQInputType_t<T> generateLLSQSystem(const std::vector<LLSQDataPoint> &data)
  {

    typename LLSQInputType_t<T>::KinematicFactors K,KWork;
    std::vector<ENSEM::EnsemComplex> vectorData;
    std::vector<LLSQDataPoint>::const_iterator it;
    bool initK = false;

    it = data.begin();

    do {
      ffKinematicFactors_t<T> genK(FormFactorDecompositionFactoryEnv::callFactory(it->matElemID));
      KWork = genK.genFactors(makeMomInvariants(it->E_f,it->E_i,it->p_f,it->p_i));

      if(it->zero.first)
	{
	  initK = true;
	  K.reDim(KWork.getB(),1,KWork.getM());
	  for(int i = 0; i < KWork.getM(); i++)
	    K.loadEnsemElement(0,i,KWork(0,i));
	  vectorData.push_back(it->zero.second);

	  if(it->one.first)
	    {
	      K.append_row(KWork.getRow(1));
	      vectorData.push_back(it->one.second);
	    }
	  if(it->two.first)
	    {
	      K.append_row(KWork.getRow(2));
	      vectorData.push_back(it->two.second);
	    }
	  if(it->three.first)
	    {
	      K.append_row(KWork.getRow(3));
	      vectorData.push_back(it->three.second);
	    }
	}
      else if(it->one.first)
	{
	  initK = true;
	  K.reDim(KWork.getB(),1,KWork.getM());
	  for(int i = 0; i < KWork.getM(); i++)
	    K.loadEnsemElement(0,i,KWork(1,i));
	  vectorData.push_back(it->zero.second);

	  if(it->two.first)
	    {
	      K.append_row(KWork.getRow(2));
	      vectorData.push_back(it->two.second);
	    }
	  if(it->three.first)
	    {
	      K.append_row(KWork.getRow(3));
	      vectorData.push_back(it->three.second);
	    }
	}
      else if(it->two.first)
	{
	  initK = true;
	  K.reDim(KWork.getB(),1,KWork.getM());
	  for(int i = 0; i < KWork.getM(); i++)
	    K.loadEnsemElement(0,i,KWork(2,i));
	  vectorData.push_back(it->zero.second);

	  if(it->three.first)
	    {
	      K.append_row(KWork.getRow(3));
	      vectorData.push_back(it->three.second);
	    }
	}
      else if(it->three.first)
	{
	  initK = true;
	  K.reDim(KWork.getB(),1,KWork.getM());
	  for(int i = 0; i < KWork.getM(); i++)
	    K.loadEnsemElement(0,i,KWork(3,i));
	  vectorData.push_back(it->zero.second);
	}
      else
	{
	  SPLASH("there was a data error in this context, basically christian is a moron");
	  exit(1);
	}
    } while(false); // execute once to force a scoped loop

    do
      {
	it ++;
	if(it == data.end())
	  break;

	ffKinematicFactors_t<T> genK(FormFactorDecompositionFactoryEnv::callFactory(it->matElemID));
	KWork = genK.genFactors(makeMomInvariants(it->E_f,it->E_i,it->p_f,it->p_i));

	if(it->zero.first)
	  {
	    K.append_row(KWork.getRow(0));
	    vectorData.push_back(it->zero.second);
	  }
	if(it->one.first)
	  {
	    K.append_row(KWork.getRow(1));
	    vectorData.push_back(it->one.second);
	  }
	if(it->two.first)
	  {
	    K.append_row(KWork.getRow(2));
	    vectorData.push_back(it->two.second);
	  }
	if(it->three.first)
	  {
	    K.append_row(KWork.getRow(3));
	    vectorData.push_back(it->three.second);
	  }

      } while(true);

    SEMBLE::SembleVector<std::complex<double> > matElems(KWork.getB(),vectorData.size());
    for(int i =0; i < vectorData.size(); i++)
      matElems.loadEnsemElement(i,vectorData[i]);

    return LLSQInputType_t<T>(K,matElems);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  template<typename T>
  struct LLSQBaseSolver_t
  {
    typedef typename ADAT::Handle<LLSQRetTypeBase_t<T> > LLSQRetTypeBase_h;
    typedef typename ADAT::Handle<LLSQInputType_t<T> > LLSQInputType_h;
    virtual LLSQRetTypeBase_h operator()(const LLSQInputType_h &) const = 0;
    virtual ~LLSQBaseSolver_t(void) {}
  };



}

#endif
