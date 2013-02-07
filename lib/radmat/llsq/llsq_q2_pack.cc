/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : llsq_q2_pack.cc

 * Purpose :

 * Creation Date : 17-01-2013

 * Last Modified : Sun Jan 27 09:26:28 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "llsq_q2_pack.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/ff/ff_gen_llsq_row.h"
#include "radmat/ff/formfactor_factory.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/splash.h"
#include "semble/semble_vector.h"
#include "semble/semble_matrix.h"
#include "semble/semble_algebra.h"
#include "semble/semble_file_management.h"
#include "adat/handle.h"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include "ensem/ensem.h"




namespace radmat
{

  namespace
  {
    void dump_zeroed(const std::vector<ENSEM::EnsemComplex> &dat, const int &idx, const LLSQDataPoint &pt, const int nbins)
    {
      std::string pth = SEMBLE::SEMBLEIO::getPath();
      std::stringstream ss;
      ss << pth << "llsq";
      SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
      ss << "/zeroed";
      SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
      ss << "/" << pt.matElemID << "_V_" << idx << "_pf" << pt.p_f[0] << pt.p_f[1] << pt.p_f[2]
        << "_pi" << pt.p_i[0] << pt.p_i[1] << pt.p_i[2] << "___zeroed.jack"; 

      ENSEM::EnsemVectorComplex e;
      e.resize(nbins);
      e.resizeObs(dat.size());

      for(unsigned int i = 0; i < dat.size(); ++i)
        ENSEM::pokeObs(e,dat[i],int(i));

      ENSEM::write(ss.str(),e);
    }

    void dump_non_zeroed(const std::vector<ENSEM::EnsemComplex> &dat, const int &idx, const LLSQDataPoint &pt, const int nbins)
    {
      std::string pth = SEMBLE::SEMBLEIO::getPath();
      std::stringstream ss;
      ss << pth << "llsq";
      SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
      ss << "/non_zeroed";
      SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
      ss << "/" << pt.matElemID << "_V_" << idx << "_pf" << pt.p_f[0] << pt.p_f[1] << pt.p_f[2]
        << "_pi" << pt.p_i[0] << pt.p_i[1] << pt.p_i[2] << "___non_zeroed.jack"; 

      ENSEM::EnsemVectorComplex e;
      e.resize(nbins);
      e.resizeObs(dat.size());

      for(unsigned int i = 0; i < dat.size(); ++i)
        ENSEM::pokeObs(e,dat[i],int(i));

      ENSEM::write(ss.str(),e);
    }

  bool isActive(const LLSQDataPoint &d)
  { 
    return !!!(!d.zero.first && !d.one.first && !d.two.first && !d.three.first);   // true if active so see if we can accumulate a false and then if we do it had active
  }



  } // anonomyous



  void LLSQDataPointQ2Pack::zeroFilter(void)
  {
    typedef std::complex<double> T;

    LLSQDataPointQ2Pack::iterator it; 
    const unsigned int dim = m_map.begin()->second.size();

    // sanity
    for(it = m_map.begin(); it != m_map.end(); ++it)
      POW2_ASSERT(dim == it->second.size()); 

    for(unsigned int pack = 0; pack < dim; ++pack)
    {
      LLSQDataPoint point(m_map.begin()->second[pack]);    
      LLSQInputType_t<T>::KinematicFactors K;
      ffKinematicFactors_t<T> genK(FormFactorDecompositionFactoryEnv::callFactory(point.matElemID));
      K = genK.genFactors(makeMomInvariants(point.E_f,point.E_i,point.p_f,point.p_i,point.mom_fac));
      SEMBLE::SembleVector<T> row; 

      row = K.getRow(0);
      row.zeros();

      bool zero,one,two,three; 

      zero = row == K.getRow(0);
      one = row == K.getRow(1);
      two = row == K.getRow(2);
      three = row == K.getRow(3); 

      std::vector<ENSEM::EnsemComplex> ez,eo,etw,eth; 

      for(it = m_map.begin(); it != m_map.end(); ++it)
      {

        // need to consider if the data was initially active -- cant just reset
        // the bool since we may not have generated the appropriate 3pt correlators
        // to be able to calculate this component
        if(it->second[pack].zero.first)
          it->second[pack].zero.first = !zero;
        if(it->second[pack].one.first)
          it->second[pack].one.first = !one;
        if(it->second[pack].two.first)
          it->second[pack].two.first = !two;
        if(it->second[pack].three.first)
          it->second[pack].three.first = !three; 

        ez.push_back(it->second[pack].zero.second);
        eo.push_back(it->second[pack].one.second);
        etw.push_back(it->second[pack].two.second);
        eth.push_back(it->second[pack].three.second);
      }

      if(zero)
        dump_zeroed(ez,0,point,K.getB());
      else
        dump_non_zeroed(ez,0,point,K.getB());

      if(one)
        dump_zeroed(eo,1,point,K.getB());
      else
        dump_non_zeroed(eo,1,point,K.getB()); 

      if(two)
        dump_zeroed(etw,2,point,K.getB());
      else
        dump_non_zeroed(etw,2,point,K.getB()); 

      if(three)
        dump_zeroed(eth,3,point,K.getB());
      else
        dump_non_zeroed(eth,3,point,K.getB()); 

    }

  }


  bool LLSQDataPointQ2Pack::haveData(void)
  {
    std::vector<LLSQDataPoint> chk = begin()->second;
    std::vector<LLSQDataPoint>::const_iterator it;
 
    for(it = chk.begin(); it != chk.end(); ++it)
      if( isActive(*it) )
        return true;
  
    return false; 
  }


}

