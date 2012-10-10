#ifndef BUILD_Q2_PACKS_H_H_GUARD
#define BUILD_Q2_PACKS_H_H_GUARD


#include "three_point.h"
#include "radmat/llsq/llsq_gen_system.h"
#include "radmat/llsq/llsq_q2_pack.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/splash.h"
#include "adat/handle.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <map>

#define  USE_OMP_LOAD_DATA_BUILDQ2PACKS

#ifdef USE_OMP_LOAD_DATA_BUILDQ2PACKS
#include <omp.h>
#endif


// #define DEBUG_NASTY_NESTED_STRUCT_WITH_POINTERS_THAT_I_HATE
// actually its not that bad

// #define DEBUGGIN_BUILD_PACKS


// parallel region in here for normalizing the correlators and such,
// the vectors already have their storage allocated, we will be 
// writing one thread per elem so its safe


namespace radmat
{


  template<typename T>
    struct BuildQ2Packs
    {
      typedef typename std::vector<ThreePointCorrelator<T> > vCor;

      BuildQ2Packs(void);


      void load(const vCor &cor);
      void normalizeZ(void);
      void normalizeExp(void);


      std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > getQ2Packs(void);


      private:
      struct pMinus
      {

        pMinus(const int xx, const int yy, const int zz)
          : x(xx) , y(yy) , z(zz)
        {}

        pMinus(const pProp_t &moms)
          : x(moms.momSource[0] - moms.momSink[0]),
          y(moms.momSource[1] - moms.momSink[1]), 
          z(moms.momSource[2] - moms.momSink[2])
        {} 


        bool equal(const pMinus &o) const
        {
          return ((x == o.x) && (y == o.y) && (z == o.z));
        }

        private:
        pMinus(void);

        public:
        int x;
        int y;
        int z;
      };

      private:
      struct SortThreePoint
      {
        typedef ThreePointCorrelator<T> C3_t;


        SortThreePoint(void)
          : zero(NULL) , one(NULL) , two(NULL) , three(NULL) 
        {  }

        SortThreePoint(const C3_t *insertion)
          : zero(NULL) , one(NULL) , two(NULL) , three(NULL)
        {
          POW2_ASSERT( insert(insertion) );
        }

        // shallow copy
        SortThreePoint(const SortThreePoint &o)
          : zero(o.zero) , one(o.one) , two(o.two) , three(o.three)
        { }

        // shallow copy
        SortThreePoint& operator=(const SortThreePoint &o)
        {
          if(this != &o)
          {
            zero = o.zero;
            one = o.one;
            two = o.two;
            three = o.three;
          }
          return *this;
        }


        bool insert(const C3_t *insertion)
        {

          if(insertion->lorentz == 0)
            if(!!!zero)
              zero = insertion;
            else
              POW2_ASSERT(false);
          else if (insertion->lorentz == 1)
            if(!!!one)
              one = insertion;
            else
              POW2_ASSERT(false);
          else if(insertion->lorentz == 2)
            if(!!!two)
              two = insertion;
            else
              POW2_ASSERT(false);
          else if(insertion->lorentz == 3)
            if(!!!three)
              three = insertion;
            else
              POW2_ASSERT(false);
          else
          {
            SPLASH("an erro occured in this context, exiting");
            POW2_ASSERT(false);
          }

          return true;
        } 

        bool insert_ref(const C3_t * insertion)
        {
          const C3_t * bar;
          if(zero) 
            bar = zero;
          else if(one)
            bar = one;
          else if(two)
            bar = two;
          else if(three) 
            bar = three;
          else
            return insert(insertion);

          if ( (bar->hel_source == insertion->hel_source) 
              && (bar->hel_sink == insertion->hel_sink) 
              && equals(bar->E_source , insertion->E_source)
              && equals(bar->E_sink , insertion->E_sink) 
              && (bar->mom.momSource == insertion->mom.momSource)
              && (bar->mom.momSink == insertion->mom.momSink))
            return insert(insertion);


          return false;
        }

        // adat wants to give back something like an obs vector of bools which is 
        // annoying so we are just going to hack an equals function here
        bool equals(const ENSEM::EnsemReal &a, const ENSEM::EnsemReal &b)
        {
          if(a.size() != b.size())
            return false;

          // even here the ensem stuff wants to be a pain.. sigh
          for(int cfg = 0; cfg < a.size(); ++cfg)
            if(SEMBLE::toScalar(a.elem(cfg)) != SEMBLE::toScalar(b.elem(cfg)))
              return false;

          return true;
        }


        const C3_t *zero;
        const C3_t *one;
        const C3_t *two;
        const C3_t *three;

      };


      private:
      void normZ(const int i);
      void normE(const int i);

      bool haveCorr;
      bool normalizedZ;
      bool normalizedExp;
      vCor m_cor;
    };


  //
  // Impl
  //////////////

  template<typename T>
    BuildQ2Packs<T>::BuildQ2Packs(void)
    : haveCorr(false) , normalizedZ(false) , normalizedExp(false)
    {}


  template<typename T>
    void BuildQ2Packs<T>::load(const vCor &cor) 
    {
      m_cor = cor;
      haveCorr = true;
      normalizedZ = false;
      normalizedExp = false;
    }


  template<typename T>
    void BuildQ2Packs<T>::normalizeZ(void)
    {
      POW2_ASSERT(haveCorr);

      std::cout << __func__ << " Normalizing.." << std::endl;

      if(!!!normalizedZ)
      {

        int i;
        int sz = m_cor.size();

#ifdef DEBUGGIN_BUILD_PACKS
        std::cout << "pre norm Z" << std::endl;
      std::cout << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(m_cor.begin()->C3pt,0))) << " +/- " 
<< std::cout << SEMBLE::toScalar(ENSEM::variance(ENSEM::peekObs(m_cor.begin()->C3pt,0))) << std::endl;
        std::cout << "Z_source" << std::endl;
        std::cout << SEMBLE::toScalar(ENSEM::mean(m_cor.begin()->Z_source)) << " +/- " 
<< SEMBLE::toScalar(ENSEM::variance(m_cor.begin()->Z_source)) << std::endl;
        std::cout << "Z_sink" << std::endl;
        std::cout << SEMBLE::toScalar(ENSEM::mean(m_cor.begin()->Z_sink)) << " +/- " 
<<  SEMBLE::toScalar(ENSEM::variance(m_cor.begin()->Z_sink)) << std::endl;
#endif

#ifdef USE_OMP_LOAD_DATA_BUILDQ2PACKS
#pragma omp parallel for shared(i,sz)
#endif
        for(i = 0; i < sz; i++)        
          normZ(i);

        /*
           std::cout << "post norm Z" << std::endl;
           std::cout << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(m_cor.begin()->C3pt,0))) << std::endl;
         */

      }
      normalizedZ = true;
    }

  // NB the 1/2E energy factor is assumed to be in the Z's
  template<typename T>
    void BuildQ2Packs<T>::normZ(const int i)
    {
      typename SEMBLE::PromoteEnsemVec<T>::Type foo = m_cor.at(i).C3pt;
      m_cor.at(i).C3pt = ( (foo)/(ENSEM::conj(m_cor.at(i).Z_source)*(m_cor.at(i).Z_sink) ) );
    }


  template<typename T>
    void BuildQ2Packs<T>::normalizeExp(void)
    {
      POW2_ASSERT(haveCorr);

      std::cout << __func__ << " removing principal time dependence" << std::endl;

      if(!!!normalizedExp)
      {

        int i;
        int sz = m_cor.size();

#ifdef DEBUGGIN_BUILD_PACKS

        std::cout << "pre norm E" << std::endl;
        std::cout << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(m_cor.begin()->C3pt,0))) << std::endl;
        std::cout << "E_source " << SEMBLE::toScalar(ENSEM::mean(m_cor.begin()->E_source)) << " +/- " 
          << SEMBLE::toScalar(ENSEM::variance(m_cor.begin()->E_source))  << std::endl;
        std::cout << "E_sink " <<  SEMBLE::toScalar(ENSEM::mean(m_cor.begin()->E_sink)) << " +/- " 
          << SEMBLE::toScalar(ENSEM::variance(m_cor.begin()->E_sink))  << std::endl;


#endif

#ifdef USE_OMP_LOAD_DATA_BUILDQ2PACKS
#pragma omp parallel for shared(i,sz)
#endif
        for(i = 0; i < sz; i++)
          normE(i);

#ifdef DEBUGGIN_BUILD_PACKS

        std::cout << "post norm E" << std::endl;
        std::cout << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(m_cor.begin()->C3pt,0))) << std::endl;
#endif

      }

      normalizedExp = true;
    }


  template<typename T>
    void BuildQ2Packs<T>::normE(const int i)
    {
      typename SEMBLE::PromoteEnsemVec<T>::Type foo = m_cor.at(i).C3pt;
      ENSEM::EnsemReal E_source = m_cor.at(i).E_source;
      ENSEM::EnsemReal E_sink = m_cor.at(i).E_sink;
      const int t_source = m_cor.at(i).t_source;
      const int t_sink = m_cor.at(i).t_sink;
      const int Lt = abs(t_source - t_sink);

      for(int t = 0; t <= Lt; t++)
      {

        ENSEM::EnsemReal factor;
        factor =  ENSEM::exp(E_sink * SEMBLE::toScalar(double(t_sink - t)))
          * ENSEM::exp(E_source * SEMBLE::toScalar(double(t - t_source) ) );
        pokeObs(m_cor.at(i).C3pt,peekObs(foo,t)*factor,t);
      }
    }


  // NB: only one type of matrix elem is assumed here, there may 
  // be problems with different irrep Q2 pieces not quite being compatible.. 
  // we will wait and see when we get real data and just redefine this function

  // we've also assumed the dispersion relation is exact to first order in p^2 there
  // could be incompatible Q2 problems here too I think..
  template<typename T>
    std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > BuildQ2Packs<T>::getQ2Packs(void)
    {
      POW2_ASSERT(normalizedZ && normalizedExp);

      std::vector<std::pair<pMinus,std::vector<int> > > sortQ2; // hold pMinus -> vector by index

      for(unsigned int i = 0; i < m_cor.size(); i++)
      {
        typename std::vector<std::pair<pMinus,std::vector<int> > >::iterator it;
        for(it = sortQ2.begin(); it != sortQ2.end(); it++)
          if(it->first.equal(pMinus(m_cor[i].mom)))
            break;

        if(it != sortQ2.end())
          it->second.push_back(i);
        else
          sortQ2.push_back(std::pair<pMinus,std::vector<int> >
              (pMinus(m_cor[i].mom),std::vector<int>(1,i)));
      }

      // now they are sorted by p_minus_spatial which effectively sorts by 
      // Q2 w/o dealing with the nasty business of defining equality 
      // between ensembles

      int packnum;
      int npacks = sortQ2.size();

#ifdef  DEBUG_NASTY_NESTED_STRUCT_WITH_POINTERS_THAT_I_HATE

      SPLASH("debuggin");
      for(unsigned int i = 0; i < sortQ2.size(); ++i)
      {
        std::cout << "pminus < " << sortQ2[i].first.x << " " 
          << sortQ2[i].first.y << " " 
          << sortQ2[i].first.z << " " << std::endl;

        for(unsigned int j = 0; j < sortQ2[i].second.size(); ++j)
          std::cout << sortQ2[i].second[j] << " " ;

        std::cout << std::endl;

      }

      SPLASH("debug");
      for(unsigned int i = 0; i < m_cor.size(); i++)
        for(int t = 0; t < m_cor[i].C3pt.numElem(); t++)  
        {
          std::cout << i  << " " << t << " " << m_cor[i].C3pt.size() 
            << " " << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(m_cor[i].C3pt,t)))
            << std::endl;

        }
#endif


      // set up the return vector and pre allocate the storage so we can go 
      // into a parallel region safely
      std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > q2_pack;  
      for(int i = 0; i < npacks; i++)
      {
        ADAT::Handle<LLSQDataPointQ2Pack> bar(new LLSQDataPointQ2Pack());
        q2_pack.push_back(bar);
      }

#ifdef USE_OMP_LOAD_DATA_BUILDQ2PACKS
#pragma omp parallel for shared(packnum,npacks)
#endif

      // loop on the sets of same Q2 and organize them 
      for(packnum = 0; packnum < npacks; packnum++)
      {

        std::vector<int> m_threepoints = sortQ2[packnum].second;
        ADAT::Handle<LLSQDataPointQ2Pack> m_handle = q2_pack[packnum];
        POW2_ASSERT(&*m_handle); // check allocation

        m_handle->setQ2(m_cor[m_threepoints[0]].Q2);

        std::vector<int>::const_iterator three_point_iterator; 
        std::vector<SortThreePoint> lorentz_packs;
        typename std::vector<SortThreePoint>::iterator lorentz_it;

        // brute force sort them into sets of 4 according to lorentz index 
        for(three_point_iterator = m_threepoints.begin();
            three_point_iterator != m_threepoints.end();
            ++three_point_iterator)
        {
          for(lorentz_it = lorentz_packs.begin();
              lorentz_it != lorentz_packs.end();
              ++lorentz_it)
            if(lorentz_it->insert_ref( &m_cor[*three_point_iterator] ) )
              break;

          if(lorentz_it == lorentz_packs.end())
            lorentz_packs.push_back(SortThreePoint(&m_cor[*three_point_iterator]));
        }

#ifdef  DEBUG_NASTY_NESTED_STRUCT_WITH_POINTERS_THAT_I_HATE

#pragma omp critical
        {
          SPLASH("debug 2");

          std::cout << "Thread " << omp_get_thread_num() << std::endl;

          std::cout << lorentz_packs.size() << std::endl;

          for(lorentz_it = lorentz_packs.begin();
              lorentz_it != lorentz_packs.end();
              ++lorentz_it)
          {
            if(!!!lorentz_it->zero)
            {
              std::cout << "caught null pointer " << std::endl;
              std::cout  << " &*(lorentz_it->0) " << &(void*)(lorentz_it->zero) << std::endl;
            }
            if(!!!lorentz_it->one)
            {
              std::cout << "caught null pointer " << std::endl;
              std::cout  << " &*(lorentz_it->1) " << &(void*)(lorentz_it->one) << std::endl;
            }
            if(!!!lorentz_it->two)
            {
              std::cout << "caught null pointer " << std::endl;
              std::cout  << " &*(lorentz_it->2) " << &(void*)(lorentz_it->two) << std::endl;
            }
            if(!!!lorentz_it->three)
            {
              std::cout << "caught null pointer " << std::endl;
              std::cout  << " &*(lorentz_it->3) " << &(void*)(lorentz_it->three) << std::endl;
            }


            for(int t = 0; t < lorentz_it->zero->C3pt.numElem(); ++t)
              std::cout << 0 << " " << t << " " 
                << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(lorentz_it->zero->C3pt,t))) << std::endl;


            for(int t = 0; t < lorentz_it->one->C3pt.numElem(); ++t)
              std::cout << 1 << " " << t << " " 
                << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(lorentz_it->one->C3pt,t))) << std::endl;

            for(int t = 0; t < lorentz_it->two->C3pt.numElem(); ++t)
              std::cout << 2 << " " << t << " " 
                << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(lorentz_it->two->C3pt,t))) << std::endl;

            for(int t = 0; t < lorentz_it->three->C3pt.numElem(); ++t)
              std::cout << 3 << " " << t << " " 
                << SEMBLE::toScalar(ENSEM::mean(ENSEM::peekObs(lorentz_it->three->C3pt,t))) << std::endl;
          }
        }
#endif


        // now  lorentz_packs is sorted to the form we need to stick into 
        // q2_pack[i] / m_handle

        int Lt = m_cor[0].C3pt.numElem();

        for(int t= 0; t < Lt; t++)
        {
          std::vector<LLSQDataPoint> llsq_data_on_t;

          for(lorentz_it = lorentz_packs.begin(); 
              lorentz_it != lorentz_packs.end();
              ++lorentz_it)
          {
            std::stringstream ss;
            LLSQDataPoint data;
            const ThreePointCorrelator<T> *writePtr;
            writePtr = NULL;

            if(lorentz_it->zero)
              writePtr = lorentz_it->zero;
            else if(lorentz_it->one)
              writePtr = lorentz_it->one;
            else if(lorentz_it->two)
              writePtr = lorentz_it->two;
            else if(lorentz_it->three)
              writePtr = lorentz_it->three;

            POW2_ASSERT(writePtr);
            ss << writePtr->elemIDBase << "_" << writePtr->hel_source 
              << "_" << writePtr->hel_sink;

            data.matElemID = ss.str();
            data.p_f = writePtr->mom.momSink;
            data.p_i = writePtr->mom.momSource;
            data.E_f = writePtr->E_sink;
            data.E_i = writePtr->E_source;
            data.mom_fac = writePtr->mom_factor;

            if(lorentz_it->zero)
            {
              data.zero.first = true;
              data.zero.second = ENSEM::peekObs(lorentz_it->zero->C3pt,t);
            }
            else
              data.zero.first = false;

            if(lorentz_it->one)
            {
              data.one.first = true;
              data.one.second = ENSEM::peekObs(lorentz_it->one->C3pt,t);
            }
            else
              data.one.first = false;  

            if(lorentz_it->two)
            {
              data.two.first = true;
              data.two.second = ENSEM::peekObs(lorentz_it->two->C3pt,t);
            }
            else
              data.two.first = false;    

            if(lorentz_it->three)
            {
              data.three.first = true;
              data.three.second = ENSEM::peekObs(lorentz_it->three->C3pt,t);
            }
            else
              data.three.first = false;


#ifdef  DEBUG_NASTY_NESTED_STRUCT_WITH_POINTERS_THAT_I_HATE
#pragma omp critical
            {
              std::cout << writePtr->mom << std::endl;
              std::cout << data.zero.first << " " << SEMBLE::toScalar(ENSEM::mean(data.zero.second)) << std::endl;
              std::cout << data.one.first << " " << SEMBLE::toScalar(ENSEM::mean(data.one.second)) << std::endl;
              std::cout << data.two.first << " " << SEMBLE::toScalar(ENSEM::mean(data.two.second)) << std::endl;
              std::cout << data.three.first << " " << SEMBLE::toScalar(ENSEM::mean(data.three.second)) << std::endl;
            }
#endif

            llsq_data_on_t.push_back(data);

          } // loop over lorentz_packs

          LLSQDataPointQ2Pack::LLSQDataPointmap_t::value_type tmp(t,llsq_data_on_t);
          m_handle->insert(tmp);

        } // loop over time


      }  // end parallel


      return q2_pack;
    }



} // namespace radmat 


#undef USE_OMP_LOAD_DATA_BUILDQ2PACKS
#undef DEBUG_NASTY_NESTED_STRUCT_WITH_POINTERS_THAT_I_HATE
#undef DEBUGGIN_BUILD_PACKS


#endif
