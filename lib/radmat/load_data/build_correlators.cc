/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : build_correlators.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Mon Dec 10 09:42:29 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "build_correlators.h"
#include "simple_world.h"
#include "generate_redstar_xml.h"  
#include "radmat_database_interface.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/llsq/llsq_gen_system.h"
#include "semble/semble_meta.h"
#include "semble/semble_key_val_db.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include <sstream>
#include <string>
#include <vector>

namespace radmat
{

  namespace
  {
    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
      {
        if(ptop.count(path) > 0)
          read(ptop,path,place);
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
            << ", path was empty, exiting" << std::endl;
          exit(1);
        }
      }

  } // namespace anonomyous 

  //! write it to a string
  std::string toString(const ThreePointCorrXMLIni_t &o)
  { 
    std::stringstream ss;
    ss << "continuumMatElemXML = " <<  o.continuumMatElemXML << "\nsource_id " 
      << o.source_id << " sink_id " << o.sink_id << " isDiagonal = " << o.isDiagonal;
    return ss.str();
  }

  //! stream it
  std::ostream& operator<<(std::ostream &o, const ThreePointCorrXMLIni_t &p)
  {
    o << toString(p);
    return o;
  }

  //! xml reader
  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrXMLIni_t &prop)
  {
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"continuumMatElemXML",prop.continuumMatElemXML,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"source_id",prop.source_id,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"sink_id",prop.sink_id,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"isDiagonal",prop.isDiagonal,__PRETTY_FUNCTION__);
  }

  //! xml writer
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrXMLIni_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"continuumMatElemXML",prop.continuumMatElemXML);
    write(xml,"source_id",prop.source_id);
    write(xml,"sink_id",prop.sink_id);
    write(xml,"isDiagonal",prop.isDiagonal);
    ADATXML::pop(xml);
  }

  std::string toString(const ThreePointCorrIni_t &prop)
  {
    std::stringstream ss;
    ss << "threePointCorrXMLIni = " << prop.threePointCorrXMLIni
      << "\nradmatDBProp = "  << prop.radmatDBProp
      << "\nmatElemID = " << prop.matElemID
      << " xi = " << prop.xi << " L_s = " << prop.L_s;  
    return ss.str();
  }

  std::ostream& operator<<(std::ostream &o, const ThreePointCorrIni_t &prop)
  {
    o << toString(prop);
    return o;
  }

  void read(ADATXML::XMLReader &xml, const std::string &path, ThreePointCorrIni_t &prop)
  { 
    ADATXML::XMLReader ptop(xml,path);
    doXMLRead(ptop,"threePointCorrXMLIni",prop.threePointCorrXMLIni,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"radmatDBProp",prop.radmatDBProp,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"matElemID",prop.matElemID,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"xi",prop.xi,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"L_s",prop.L_s,__PRETTY_FUNCTION__);
  }

  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrIni_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"threePointCorrXMLIni",prop.threePointCorrXMLIni);
    write(xml,"radmatDBProp",prop.radmatDBProp);
    write(xml,"matElemID",prop.matElemID); 
    write(xml,"xi",prop.xi);
    write(xml,"L_s",prop.L_s); 
    ADATXML::pop(xml);
  }

} // namespace radmat



// BUILD CORRELATORS
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

namespace radmat
{



  namespace
  {

    typedef radmatAllConfDatabaseInterface< Hadron::KeyHadronNPartNPtCorr_t,
            ENSEM::EnsemVectorComplex,
            SEMBLE::SembleExtendedKeyHadronNPartIrrep_t,
            SEMBLE::SembleMassOverlapData_t> DatabaseInterface_t;


    // hurray a continuum three point correlator..
    struct accumulator
    {

      typedef redstarCartMatElem::redstarCartMatElemLorentzComponent data_t;

      accumulator(void)
        : active(false) , success(true) , momentum_factor(0.)
      {
        ENSEM::EnsemComplex Ezero;
        Ezero.resize(1);
        Ezero = SEMBLE::toScalar(std::complex<double>(0.,0.));
        corr.resize(1);
        corr.resizeObs(1);
        corr = SEMBLE::toScalar(std::complex<double>(0.,0.));
        E_f = ENSEM::real(Ezero);
        E_i = E_f; 

      }

      void eat(const data_t &data, const ThreePointCorrIni_t &ini)
      {
        // break early if not active
        active = data.is_active();
        if(!!!active)
          return;

        m_ini = ini; 

        DatabaseInterface_t db(m_ini.radmatDBProp);

        // accumulate a bad list
        data_t::const_iterator it; 
        for(it = data.begin(); it != data.end(); ++it)
        {
          if(!!!db.exists(it->m_obj.redstar_xml))
            m_bad_corrs.push_back(it->m_obj.redstar_xml);
          if(!!!db.exists(it->m_obj.source_normalization))
            m_bad_norms.push_back(it->m_obj.source_normalization);
          if(!!!db.exists(it->m_obj.sink_normalization))
            m_bad_norms.push_back(it->m_obj.sink_normalization);
        }

        // break early if not successful 
        if((!!!m_bad_corrs.empty()) && (m_bad_norms.empty()))
          success = false; 

        if(!!!success)
          return;

        // set up a zero valued three point correlator
        it = data.begin();
        corr = db.fetch(it->m_obj.redstar_xml);
        corr = SEMBLE::toScalar(0.); 
        ThreePtPropagationFactor<double> propagation_factor;

        // initialize some variables
        momentum_factor = 1./ini.xi * 2.*acos(-1.)/ini.L_s;      
        p_f = it->m_obj.redstar_xml.npoint[3].irrep.mom; 
        p_i = it->m_obj.redstar_xml.npoint[1].irrep.mom;

        SEMBLE::SembleMassOverlapData_t source_tmp = db.fetch(it->m_obj.source_normalization); 

        E_f = source_tmp.E();
        E_f = SEMBLE::toScalar(0.);
        E_i = E_f;
        int ct = 0;

        for(it = data.begin(); it != data.end(); ++it)
        {
          ENSEM::EnsemVectorComplex corr_tmp = db.fetch(it->m_obj.redstar_xml);
          SEMBLE::SembleMassOverlapData_t source = db.fetch(it->m_obj.source_normalization); 
          SEMBLE::SembleMassOverlapData_t sink = db.fetch(it->m_obj.sink_normalization); 

          // the hadron key uses 1 based arrays
          const int t_source(it->m_obj.redstar_xml.npoint[1].t_slice);
          const int t_sink(it->m_obj.redstar_xml.npoint[3].t_slice); 

          for(int t_ins = 0; t_ins < corr_tmp.size(); ++t_ins)
          {
            ENSEM::EnsemReal prop = propagation_factor(sink.E(),sink.Z(),t_sink,t_ins,
                source.E(),source.Z(),t_source);
            ENSEM::pokeObs(corr_tmp,ENSEM::peekObs(corr_tmp,t_ins)/prop,t_ins);
          }

          E_f = E_f + sink.E();
          E_i = E_i + source.E();

          corr = corr + it->m_coeff*corr_tmp; 

          ++ct;
        }

        // fill in the rest of the relevant data
        E_f = E_f/SEMBLE::toScalar(double(ct));
        E_i = E_i/SEMBLE::toScalar(double(ct));
      }

      bool active; 
      bool success; 
      ThreePointCorrIni_t m_ini;
      ENSEM::EnsemVectorComplex corr; 
      ENSEM::EnsemReal E_f;
      ENSEM::EnsemReal E_i; 
      ADATXML::Array<int> p_f;
      ADATXML::Array<int> p_i; 
      double momentum_factor; // 1/xi * 2pi/L_s -- the "unit" size  
      std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
      std::vector<SEMBLE::SembleExtendedKeyHadronNPartIrrep_t> m_bad_norms; 
    };



    struct CartesianMatrixElement
    {
      CartesianMatrixElement(void) {}
      CartesianMatrixElement(const ThreePointCorrIni_t &ini, const simpleWorld::ContinuumMatElem &elem)
        : have_active_data(false) , success(true) , m_ini(ini) , m_elem(elem)
      {
        redstarCartMatElem celem(elem,ini.threePointCorrXMLIni.source_id,ini.threePointCorrXMLIni.sink_id);
        t.eat(celem.get_component(0),m_ini);
        x.eat(celem.get_component(1),m_ini);
        y.eat(celem.get_component(2),m_ini);
        z.eat(celem.get_component(3),m_ini);

        check(); 

        std::stringstream ss; 
        ss << m_ini.matElemID << "_" << m_elem.source.state.H << "_" << m_elem.sink.state.H;
        m_mat_elem_id = ss.str(); 

      }

      void check(void)
      {
        if (t.active && t.success)
        {
          E_f = t.E_f; 
          E_i = t.E_i;
          p_f = t.p_f;
          p_i = t.p_i;
          have_active_data = true; 
          tmax = t.corr.numElem(); 
        }      
        else  if (x.active && x.success)
        {
          E_f = x.E_f; 
          E_i = x.E_i;
          p_f = x.p_f;
          p_i = x.p_i;
          have_active_data = true; 
          tmax = x.corr.numElem();
        }       
        else  if (y.active && y.success)
        {
          E_f = y.E_f; 
          E_i = y.E_i;
          p_f = y.p_f;
          p_i = y.p_i;
          have_active_data = true; 
          tmax = y.corr.numElem();
        }       
        else  if (z.active && z.success)
        {
          E_f = t.E_f; 
          E_i = t.E_i;
          p_f = t.p_f;
          p_i = t.p_i;
          have_active_data = true;
          tmax = z.corr.numElem();  
        }

        if(!!!(t.success && x.success && y.success && z.success))
        {
          m_bad_corrs.reserve(t.m_bad_corrs.size() + x.m_bad_corrs.size() + y.m_bad_corrs.size() + z.m_bad_corrs.size());
          m_bad_corrs.insert(m_bad_corrs.end(),t.m_bad_corrs.begin(),t.m_bad_corrs.end());
          m_bad_corrs.insert(m_bad_corrs.end(),x.m_bad_corrs.begin(),x.m_bad_corrs.end());
          m_bad_corrs.insert(m_bad_corrs.end(),y.m_bad_corrs.begin(),y.m_bad_corrs.end());
          m_bad_corrs.insert(m_bad_corrs.end(),z.m_bad_corrs.begin(),z.m_bad_corrs.end());

          m_bad_norms.reserve(t.m_bad_norms.size() + x.m_bad_norms.size() + y.m_bad_norms.size() + z.m_bad_norms.size());
          m_bad_norms.insert(m_bad_norms.end(),t.m_bad_norms.begin(),t.m_bad_norms.end());
          m_bad_norms.insert(m_bad_norms.end(),x.m_bad_norms.begin(),x.m_bad_norms.end());
          m_bad_norms.insert(m_bad_norms.end(),y.m_bad_norms.begin(),y.m_bad_norms.end());
          m_bad_norms.insert(m_bad_norms.end(),z.m_bad_norms.begin(),z.m_bad_norms.end());

          success = false;
        }
      }

      LLSQDataPoint getLLSQPoint(const int t_ins) const
      {
        LLSQDataPoint ret; 
        ret.matElemID = m_mat_elem_id; 
        ret.p_f = p_f;
        ret.p_i = p_i;
        ret.E_f = E_f;
        ret.E_i = E_i; 
        ret.mom_fac = 1./m_ini.xi * 2.*acos(-1.)/m_ini.L_s;     
        ENSEM::EnsemComplex EZERO;
        EZERO.resize(E_f.size());
        EZERO = SEMBLE::toScalar(std::complex<double>(0.,0.));
        std::pair<bool,ENSEM::EnsemComplex> badpt(false,EZERO);

        if(t.active && t.success)
          ret.zero = std::pair<bool,ENSEM::EnsemComplex>(true,ENSEM::peekObs(t.corr,t_ins));
        else
          ret.zero = badpt;

        if(x.active && x.success)
          ret.one = std::pair<bool,ENSEM::EnsemComplex>(true,ENSEM::peekObs(x.corr,t_ins));
        else
          ret.one = badpt;

        if(y.active && y.success)
          ret.two = std::pair<bool,ENSEM::EnsemComplex>(true,ENSEM::peekObs(y.corr,t_ins));
        else
          ret.two = badpt;

        if(z.active && z.success)
          ret.three = std::pair<bool,ENSEM::EnsemComplex>(true,ENSEM::peekObs(z.corr,t_ins));
        else
          ret.three = badpt;

        return ret;      
      }


      bool have_active_data;
      bool success; 
      ENSEM::EnsemReal E_f,E_i; 
      ADATXML::Array<int> p_f, p_i; 
      ThreePointCorrIni_t m_ini; 
      simpleWorld::ContinuumMatElem m_elem; 
      std::string m_mat_elem_id; 
      accumulator t,x,y,z;
      int tmax;
      std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
      std::vector<SEMBLE::SembleExtendedKeyHadronNPartIrrep_t> m_bad_norms; 


    }; // cartesianMatrixElement


    std::vector<CartesianMatrixElement> 
      getCartesianMatrixElements(const std::vector<simpleWorld::ContinuumMatElem> & elems, const ThreePointCorrIni_t &ini)
      {
        std::vector<CartesianMatrixElement> ret;
        std::vector<simpleWorld::ContinuumMatElem>::const_iterator it; 
        for(it = elems.begin(); it != elems.end(); ++it)
          ret.push_back(CartesianMatrixElement(ini,*it)); 
        return ret; 
      }


    int dot(const ADATXML::Array<int> &a, const ADATXML::Array<int> &b)
    {
      return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }


    struct Q2HashKey
    {
      Q2HashKey(const int pf, const int pi, const int d) 
        : pf2(pf) , pi2(pi) , dotprod(d) {}

      Q2HashKey(const CartesianMatrixElement &elem)
      {
        pf2 = dot(elem.p_f,elem.p_f);
        pi2 = dot(elem.p_i,elem.p_i);
        dotprod = dot(elem.p_f,elem.p_i);
      }

      int pf2,pi2,dotprod;
    };


    bool operator<(const Q2HashKey &a, const Q2HashKey &b)
    {
      if (a.pf2 < b.pf2)
        return true;
      if(a.pi2 < b.pi2)
        return true;
      if(a.dotprod < b.dotprod)
        return true;
      return false;
    }

    bool operator>(const Q2HashKey &a, const Q2HashKey &b)
    {
      if (a.pf2 > b.pf2)
        return true;
      if(a.pi2 > b.pi2)
        return true;
      if(a.dotprod > b.dotprod)
        return true;
      return false;
    }

    bool operator==(const Q2HashKey &a, const Q2HashKey &b)
    {
      return  (!(a < b) && !(a > b));
    }

    struct compareQ2HashKey
    {

      compareQ2HashKey(const bool isDiagonal)
        : m_isDiagonal(isDiagonal) {}

      bool operator()(const Q2HashKey &a, const Q2HashKey &b) const 
      {
        if(m_isDiagonal)
        {
          // this deals with th symmetry when intial and final momenta are reversed 
          if((a.pf2 == b.pi2)&&(a.pi2 == b.pf2)&&(a.dotprod == b.dotprod))  
            return true;

          // this deals with Q2 = 0 being accessible multiple ways -- basically check that the momenta were colinear w/ same direction (q_space = 0)
          if((a.pf2 == a.pi2)&&(a.pf2 == a.dotprod)&&(b.pf2 == b.pi2)&&(b.pf2 == b.dotprod))
            return true;
        }

        return a == b; 
      }

      bool m_isDiagonal;
    };


    struct manageWork
    {
      manageWork(const ThreePointCorrIni_t & ini)
        : m_ini(ini) 
      {

        std::vector<CartesianMatrixElement> unsorted;
        unsorted = getCartesianMatrixElements(simpleWorld::getContinuumMatElemFromXML(m_ini.threePointCorrXMLIni.continuumMatElemXML),m_ini);
        sort(unsorted,m_ini.threePointCorrXMLIni.isDiagonal);
      }

      ~manageWork(void)
      {

        // dump the accumulated bad lists 
        if(!!!m_bad_corrs.empty())
        {
          ADATXML::XMLBufferWriter corrs;
          ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;
          bc.resize(m_bad_corrs.size());
          for(unsigned int i = 0; i < m_bad_corrs.size(); ++i)
            bc[i] = m_bad_corrs[i];
          write(corrs,"BadCorrs",bc);
          std::ofstream out("missing_three_point_correlators.xml");
          corrs.print(out);
          out.close();
        }

        if(!!!m_bad_norms.empty())
        {
          ADATXML::XMLBufferWriter norms;
          ADATXML::Array<SEMBLE::SembleExtendedKeyHadronNPartIrrep_t> bn;
          bn.resize(m_bad_norms.size());
          for(unsigned int i = 0; i < m_bad_norms.size(); ++i)
            bn[i] = m_bad_norms[i];
          write(norms,"BadNorms",bn);
          std::ofstream out("missing_normalizations.xml");
          norms.print(out);
          out.close(); 
        }
      }


      // brute force sort using momenta 
      void sort(const std::vector<CartesianMatrixElement> &unsorted, const bool isDiagonal)
      {   
        // this deals with some stupid ambiguity associated with the "Anything that *might be* a declaration *is* a declaration" idea
        // the result of which was a headache and realization that one can't construct the comparator inplace..
        compareQ2HashKey my_comparator(isDiagonal);
        std::map<Q2HashKey,std::vector<CartesianMatrixElement>, compareQ2HashKey > my_map(my_comparator);
        typedef std::map<Q2HashKey,std::vector<CartesianMatrixElement>, compareQ2HashKey >::value_type value_type; 
        typedef std::map<Q2HashKey,std::vector<CartesianMatrixElement>, compareQ2HashKey >::iterator my_map_iterator_t; 

        std::vector<CartesianMatrixElement>::const_iterator it;
        my_map_iterator_t map_it; 
        
        for(it = unsorted.begin(); it != unsorted.end(); ++it)
        {
          if(!!!it->success)
            push_bad_list(*it);
          if(!!!it->have_active_data)
            continue;

          map_it = my_map.find(Q2HashKey(*it));

          if(map_it == my_map.end())
          {
            std::vector<CartesianMatrixElement> dum;
            dum.push_back(*it);
            my_map.insert(value_type(Q2HashKey(*it),dum));
          }
          else
            map_it->second.push_back(*it);
        }

        sorted_by_Q2.clear(); 
        std::vector<CartesianMatrixElement>::const_iterator vector_iterator; 
        for(map_it = my_map.begin(); map_it != my_map.end(); ++map_it)
        {
          std::vector<CartesianMatrixElement> tmp; 

          for(vector_iterator = map_it->second.begin(); 
              vector_iterator != map_it->second.end(); 
              ++vector_iterator)
            tmp.push_back(*it); 

          sorted_by_Q2.push_back(tmp); 
        }
      }


      void push_bad_list(const CartesianMatrixElement &elem)
      {
        m_bad_corrs.reserve(m_bad_corrs.size() + elem.m_bad_corrs.size());
        m_bad_corrs.insert(m_bad_corrs.end(),elem.m_bad_corrs.begin(),elem.m_bad_corrs.end());
        m_bad_norms.reserve(m_bad_norms.size() + elem.m_bad_norms.size());
        m_bad_norms.insert(m_bad_norms.end(),elem.m_bad_norms.begin(),elem.m_bad_norms.end()); 
      }

      std::vector<std::vector<CartesianMatrixElement> >
        get_sorted_data(void) const
        {
          return sorted_by_Q2; 
        }

      // data 
      ThreePointCorrIni_t m_ini; 
      std::vector<std::vector<CartesianMatrixElement> > sorted_by_Q2; 
      std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
      std::vector<SEMBLE::SembleExtendedKeyHadronNPartIrrep_t> m_bad_norms; 

    }; // manageWork




  } // namespace anonomyous 


  std::vector<ADAT::Handle<LLSQDataPointQ2Pack> >
    BuildCorrelators::build_correlators(void)
    {
      POW2_ASSERT(have_ini);
      manageWork work_manager(m_ini);
      std::vector<std::vector<CartesianMatrixElement> > sorted_by_q2(work_manager.get_sorted_data());
      std::vector<std::vector<CartesianMatrixElement> >::const_iterator it_q2;
      std::vector<CartesianMatrixElement>::const_iterator it_pack; 
      std::vector<ADAT::Handle<LLSQDataPointQ2Pack> > ret; 

      // this is looping over all of the Q2 values
      for(it_q2 = sorted_by_q2.begin(); it_q2 != sorted_by_q2.end(); ++it_q2)
      {
        // each Q2 value gets its own pack
        ADAT::Handle<LLSQDataPointQ2Pack> tmp(new LLSQDataPointQ2Pack());


        // a pack is a map with key t_ins and data vector<LLSQDataPoint>
        const int tmax = it_q2->begin()->tmax;

        // so loop t_ins
        for(int t_ins = 0; t_ins < tmax; ++t_ins)
        {
          // at each t_ins make a vector of all the points in the linear system 
          std::vector<LLSQDataPoint> points_at_this_t_ins; 
          for(it_pack = it_q2->begin(); it_pack != it_q2->end(); ++it_pack)
            points_at_this_t_ins.push_back(it_pack->getLLSQPoint(t_ins));  

          tmp->insert(LLSQDataPointQ2Pack::value_type(t_ins,points_at_this_t_ins)); 
        }

        ret.push_back(tmp);
      }

      return ret; 
    }


}




