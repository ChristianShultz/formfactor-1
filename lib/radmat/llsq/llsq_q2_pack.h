#ifndef LLSQ_Q2_PACK_H_H_GUARD
#define LLSQ_Q2_PACK_H_H_GUARD


#include "llsq_gen_system.h"
#include "semble/semble_vector.h"
#include "semble/semble_meta.h"
#include <map>
#include "adat/handle.h"

namespace radmat
{

  struct LLSQDataPointQ2Pack
  {
    typedef std::map<int, std::vector<LLSQDataPoint> > LLSQDataPointmap_t;
    // map index is time of the inserted current, 
    //the vector is all the points at a given Q2

    typedef LLSQDataPointmap_t::iterator iterator;
    typedef LLSQDataPointmap_t::const_iterator const_iterator;
    typedef LLSQDataPointmap_t::value_type value_type;


    // stl wrappers

    iterator begin(void) {return m_map.begin();}
    const_iterator begin(void) const { return m_map.begin();}
    iterator end(void) {return m_map.end();}
    const_iterator end(void) const {return m_map.end();}
    void insert(const value_type &v) {m_map.insert(v);}


    // once again we need to initialize the ensems
    LLSQDataPointQ2Pack(void)
    {
      m_Q2.resize(1);
      m_Q2 = SEMBLE::toScalar(double(10000000.)); 
    }


    void zeroFilter(void);

    bool haveData(void);     

    ENSEM::EnsemReal Q2(void) const {return m_Q2;}
    void setQ2(const ENSEM::EnsemReal &Q2) {m_Q2 = Q2;}

    
    LLSQDataPointmap_t m_map;
    ENSEM::EnsemReal m_Q2;

  };



  template<typename T>
    struct LLSQRet_t_Q2Pack
    {
      typedef typename LLSQBaseSolver_t<T>::LLSQRetTypeBase_h LLSQRetTypeBase_h; 
      typedef typename std::map<int, SEMBLE::SembleVector<T> > map_t;
      // map time to the extracted vector of ffs

      typedef typename map_t::iterator iterator;
      typedef typename map_t::const_iterator const_iterator;
      typedef typename map_t::value_type value_type;


      LLSQRet_t_Q2Pack(void) {}


      LLSQRet_t_Q2Pack(const LLSQRet_t_Q2Pack &o)
        : m_map(o.m_map) , m_Q2(o.m_Q2)
      {}

      ~LLSQRet_t_Q2Pack(void) {}

      // stl wrappers
      iterator begin(void) {return m_map.begin();}
      const_iterator begin(void) const {return m_map.begin();}
      iterator end(void) {return m_map.end();}
      const_iterator end(void) const {return m_map.end();}
      ENSEM::EnsemReal Q2(void) const {return m_Q2;}
      void setQ2(const ENSEM::EnsemReal &Q2) {m_Q2 = Q2;}


      unsigned int size(void) const {return m_map.size();}

      void insert(const int t, const LLSQRetTypeBase_h &v)
      {

        if(m_map.find(t) != m_map.end())
        {
          SPLASH("map key already present, exiting");
          exit(1);
        }

        m_map.insert(value_type(t,v->m_FF));
      }

      map_t m_map;
      ENSEM::EnsemReal m_Q2;


    };


  template<typename T>
    struct LLSQRet_ff_Q2Pack
    {

      LLSQRet_ff_Q2Pack(void)
      {
        m_Q2.resize(1);
        m_Q2 = SEMBLE::toScalar(double(10000000));
      }

      LLSQRet_ff_Q2Pack(const LLSQRet_ff_Q2Pack &o)
        : m_map(o.m_map) , m_Q2(o.m_Q2)
      {
        m_Q2.resize(1);
        m_Q2 = SEMBLE::toScalar(double(10000000));
      }

      ~LLSQRet_ff_Q2Pack(void) {}


      typedef typename std::map<int, typename SEMBLE::PromoteEnsemVec<T>::Type > map_t;

      // map index to ensem FF(t)

      typedef typename map_t::iterator iterator;
      typedef typename map_t::const_iterator const_iterator;
      typedef typename map_t::value_type value_type;

      // stl wrappers
      iterator begin(void) {return m_map.begin();}
      const_iterator begin(void) const {return m_map.begin();}
      iterator end(void) {return m_map.end();}
      const_iterator end(void) const {return m_map.end();}
      ENSEM::EnsemReal Q2(void) const {return m_Q2;}
      void setQ2(const ENSEM::EnsemReal &Q2) {m_Q2 = Q2;}


      unsigned int size(void) const {return m_map.size();}

      void insert(const int ff_num, const typename SEMBLE::PromoteEnsemVec<T>::Type &ensem)
      {

        // no double insertions, this is a one off class
        if(m_map.find(ff_num) != m_map.end())
        {
          SPLASH("map key already present, exiting");
          exit(1);
        }

        m_map.insert(value_type(ff_num,ensem));
      }

      map_t m_map;
      ENSEM::EnsemReal m_Q2;

    };


  template<typename T>
    struct Q2Pack
    {
      LLSQDataPointQ2Pack latticeData;   
      LLSQRet_ff_Q2Pack<T> formFactors;
    };



  // transform from the LLSQ storage scheme to a schem indexed by ff#


  template<typename T>
    ADAT::Handle<LLSQRet_ff_Q2Pack<T> >  transformLLSQRetPack(const LLSQRet_t_Q2Pack<T> &in)
    {


      ADAT::Handle<LLSQRet_ff_Q2Pack<T> > out(new LLSQRet_ff_Q2Pack<T>);      
      typename SEMBLE::PromoteEnsemVec<T>::Type zero;

      int ncfg = in.begin()->second.getB();
      int nff = in.begin()->second.getN();      // there are some possibly bad assumptions
      int Lt = in.m_map.size();          // about structure of inputs here...

      zero.resize(ncfg);
      zero.resizeObs(Lt);
      zero = SEMBLE::toScalar(T(0));

      //      // init the map to zero
      //      for(int i = 0; i < nff; i++)
      //        out.insert(i,zero);

      typename LLSQRet_t_Q2Pack<T>::const_iterator it;

      for(int i = 0; i < nff; i++)
      {
        typename SEMBLE::PromoteEnsemVec<T>::Type tmp;
        tmp = zero;

        for(it = in.begin(); it != in.end(); it++)
          ENSEM::pokeObs(tmp,it->second.getEnsemElement(i),it->first);

        out->insert(i,tmp);
      }

      out->setQ2(in.Q2());

      return out;
    }



} // close radmat






#endif
