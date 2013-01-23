/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : build_correlators.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Wed Jan 23 13:06:17 2013

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
#include "semble/semble_file_management.h"
#include "radmat_overlap_key_val_db.h"
#include "hadron/ensem_filenames.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"
#include "jackFitter/plot.h"
#include "ensem/ensem.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <set>

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
      << o.source_id << " sink_id " << o.sink_id << " isDiagonal = " << o.isDiagonal
      << " isProjected = " << o.isProjected;
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
    doXMLRead(ptop,"isProjected",prop.isProjected,__PRETTY_FUNCTION__); 
  }

  //! xml writer
  void write(ADATXML::XMLWriter &xml, const std::string &path, const ThreePointCorrXMLIni_t &prop)
  {
    ADATXML::push(xml,path);
    write(xml,"continuumMatElemXML",prop.continuumMatElemXML);
    write(xml,"source_id",prop.source_id);
    write(xml,"sink_id",prop.sink_id);
    write(xml,"isDiagonal",prop.isDiagonal);
    write(xml,"isProjected",prop.isProjected);
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
            RadmatExtendedKeyHadronNPartIrrep_t,
            RadmatMassOverlapData_t> DatabaseInterface_t;


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
        if(!!!m_bad_corrs.empty())
          success = false; 

        if(!!!m_bad_norms.empty())
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

        RadmatMassOverlapData_t source_tmp = db.fetch(it->m_obj.source_normalization); 

        E_f = source_tmp.E();
        E_f = SEMBLE::toScalar(0.);
        E_i = E_f;
        int ct = 0;

        for(it = data.begin(); it != data.end(); ++it)
        {
          ENSEM::EnsemVectorComplex corr_tmp = db.fetch(it->m_obj.redstar_xml);
          RadmatMassOverlapData_t source = db.fetch(it->m_obj.source_normalization); 
          RadmatMassOverlapData_t sink = db.fetch(it->m_obj.sink_normalization); 

          // the hadron key uses 1 based arrays

          // NB: assumption that npt is organized like <sink, ins , source>

          const int t_source(it->m_obj.redstar_xml.npoint[3].t_slice);
          const int t_sink(it->m_obj.redstar_xml.npoint[1].t_slice); 

          
          POW2_ASSERT(t_source < t_sink); 


          for(int t_ins = t_source; t_ins <= t_sink; ++t_ins)
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
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> m_bad_norms; 
    };

    namespace
    {
      std::string append_momentum(const std::string &name, const ADATXML::Array<int> &mom)
      {
        std::stringstream ss; 
        ADATXML::Array<int> can_mom = FF::canonicalOrder(mom);
        ss << name << "_p" << can_mom[0] << can_mom[1] << can_mom[2];
        return ss.str();  
      }


      // adapted from something i found on stack overflow
      // pair<iterator,bool> insert (const value_type& val);
      // return a pair, with its member pair::first set to an iterator pointing to either the newly inserted 
      // element or to the equivalent element already in the set. The pair::second element in the pair is set
      // to true if a new element was inserted or false if an equivalent element already existed.
      // 
      // since we don't have a comparator we will use the string output of the keys for simplicity
      // I think this should still be unique but didn't check carefully
      template<typename T>
        struct NotDuplicate
        {
          bool operator()(const T &elem)
          {
            std::stringstream ss; 
            ss << elem;
            return s_.insert(ss.str()).second; 
          }

          private: 
          std::set<std::string> s_; 
        };

    }

    struct CartesianMatrixElement
    {
      CartesianMatrixElement(void) 
      {
        // cant copy around uninitialized ensems and std vector
        // loves to move around empty objects so fake up and ensem
        // to avoid the checkResize error in ENSEM
        E_f.resize(1);
        E_f = SEMBLE::toScalar(double(0.)); 
        E_i = E_f; 
      }


      CartesianMatrixElement(const ThreePointCorrIni_t &ini, const simpleWorld::ContinuumMatElem &elem)
        : have_active_data(false) , success(true) , m_ini(ini) , m_elem(elem)
      {


        // hardwire some assumptions about projected operator naming here.. if something changes we can just build some
        // sort of hash to retreive the correct name but the current convetions should be good for a bit
        if(m_ini.threePointCorrXMLIni.isProjected)
        {
          m_elem.source.state.name = append_momentum(m_elem.source.state.name,m_elem.source.state.mom);
          m_elem.sink.state.name = append_momentum(m_elem.sink.state.name,m_elem.sink.state.mom);
        }

        redstarCartMatElem celem(m_elem,ini.threePointCorrXMLIni.source_id,ini.threePointCorrXMLIni.sink_id);
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

#if 0
          std::cout << __PRETTY_FUNCTION__ << std::endl;
          std::cout << "size <" << t.m_bad_corrs.size() << " " << x.m_bad_corrs.size()
            << " " << y.m_bad_corrs.size() << " " << z.m_bad_corrs.size() << std::endl;

          std::vector<Hadron::KeyHadronNPartNPtCorr_t >::const_iterator corr_it;
          std::vector<RadmatExtendedKeyHadronNPartIrrep_t>::const_iterator norm_it; 
          const std::vector<Hadron::KeyHadronNPartNPtCorr_t >* corr_ptr;
          const std::vector<RadmatExtendedKeyHadronNPartIrrep_t>* norm_ptr; 
          std::string idx;

          idx = "t";
          corr_ptr = &t.m_bad_corrs; 
          if(!!!corr_ptr->empty())
          {
            std::cout << idx << std::endl;
            for(corr_it = corr_ptr->begin(); corr_it != corr_ptr->end(); ++corr_it)
              std::cout << *corr_it << std::endl;
          }

          idx = "x";
          corr_ptr = &x.m_bad_corrs; 
          if(!!!corr_ptr->empty())
          {
            std::cout << idx << std::endl;
            for(corr_it = corr_ptr->begin(); corr_it != corr_ptr->end(); ++corr_it)
              std::cout << *corr_it << std::endl;
          }


          idx = "y";
          corr_ptr = &y.m_bad_corrs; 
          if(!!!corr_ptr->empty())
          {
            std::cout << idx << std::endl;
            for(corr_it = corr_ptr->begin(); corr_it != corr_ptr->end(); ++corr_it)
              std::cout << *corr_it << std::endl;
          }


          idx = "z";
          corr_ptr = &z.m_bad_corrs; 
          if(!!!corr_ptr->empty())
          {
            std::cout << idx << std::endl;
            for(corr_it = corr_ptr->begin(); corr_it != corr_ptr->end(); ++corr_it)
              std::cout << *corr_it << std::endl;
          }


#endif


          success = false;
        }

        // if we don't have any active data we need to initialize the 
        // ensembles to avoid checkResize error in the event that we copy
        // this non active data matrix element around for some reason
        if(!!!have_active_data)
        {
          E_f = t.E_f;
          E_i = t.E_i; 
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
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> m_bad_norms; 


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


    double Mink_qsq(const CartesianMatrixElement & elem, const double E_f , const double E_i, const double factor)
    { 
      ADATXML::Array<double> p_f,p_i; 
      p_f.resize(3);
      p_i.resize(3); 
      double p_i_sq(0.);
      double p_f_sq(0.); 
      double q_space_sq(0.);

      for(int i = 0; i < 3; ++i)
      {
        p_f[i] = factor * double(elem.p_f[i]); 
        p_i[i] = factor * double(elem.p_i[i]);
        p_f_sq += p_f[i]*p_f[i];
        p_i_sq += p_i[i]*p_i[i];
        q_space_sq += (p_i[i] - p_f[i])*(p_i[i] - p_f[i]);
      }

      double q_time =  sqrt(E_i*E_i + p_i_sq) - sqrt(E_f*E_f + p_f_sq);  

      return -((q_time*q_time) - q_space_sq); 
    }


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

        write_projected_matrix_elements();
        dump_xml();
      }


      void dump_xml(void)
      {
        // dump the accumulated bad lists 
        if(!!!m_bad_corrs.empty())
        {
          ADATXML::XMLBufferWriter corrs;
          ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> bc;

          std::vector<Hadron::KeyHadronNPartNPtCorr_t> seen;
          std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it;
          NotDuplicate<Hadron::KeyHadronNPartNPtCorr_t> predicate; 

          for(it = m_bad_corrs.begin(); it != m_bad_corrs.end(); ++it)
            if (predicate(*it))
              seen.push_back(*it); 

          bc.resize(seen.size()); 
          for(unsigned int i = 0; i < seen.size(); ++i)
            bc[i] = seen[i];
          write(corrs,"NPointList",bc);

          std::string pth = SEMBLE::SEMBLEIO::getPath();
          SEMBLE::SEMBLEIO::makeDirectoryPath(pth + std::string("missing")); 

          std::ofstream out("missing/npt.list.xml");
          corrs.print(out);
          out.close();


          out.open("missing/npt.ensemFileNames.list"); 
          for(it = seen.begin(); it != seen.end(); ++it)
            out << Hadron::ensemFileName(*it) << "\n";
          out.close(); 

        }

        if(!!!m_bad_norms.empty())
        {
          ADATXML::XMLBufferWriter norms;
          ADATXML::Array<RadmatExtendedKeyHadronNPartIrrep_t> bn;

          std::vector<RadmatExtendedKeyHadronNPartIrrep_t> seen;
          std::vector<RadmatExtendedKeyHadronNPartIrrep_t>::const_iterator it;
          NotDuplicate<RadmatExtendedKeyHadronNPartIrrep_t> predicate; 

          for(it = m_bad_norms.begin(); it != m_bad_norms.end(); ++it)
            if (predicate(*it))
              seen.push_back(*it); 

          bn.resize(seen.size());
          for(unsigned int i = 0; i < seen.size(); ++i)
            bn[i] = seen[i];

          write(norms,"BadNorms",bn);

          std::string pth = SEMBLE::SEMBLEIO::getPath();
          SEMBLE::SEMBLEIO::makeDirectoryPath(pth + std::string("missing")); 

          std::ofstream out("missing/normalizations.list.xml");
          norms.print(out);
          out.close(); 
        }

      }




      void sort(const std::vector<CartesianMatrixElement> &unsorted, const bool isDiagonal)
      {

        double E_f, E_i; // rest energies..       
        std::vector<CartesianMatrixElement>::const_iterator it;
        bool found_f(false) , found_i(false); 

        for(it = unsorted.begin(); it != unsorted.end(); ++it)
        {

          if(!!!it->success)
            push_bad_list(*it);

          if(!!!it->have_active_data)
            continue;

          if(!!!found_f)
          {
            ADATXML::Array<int> p_f = it->p_f;
            if( (p_f[0] == 0) && (p_f[1] == 0) && (p_f[2] == 0) )
            {
              E_f = SEMBLE::toScalar(ENSEM::mean(it->E_f)); 
              found_f = true; 
            }
          }


          if(!!!found_i)
          {
            ADATXML::Array<int> p_i = it->p_i;
            if( (p_i[0] == 0) && (p_i[1] == 0) && (p_i[2] == 0) )
            {
              E_i = SEMBLE::toScalar(ENSEM::mean(it->E_i)); 
              found_i = true; 
            }

          }

        } // for it


        if(!!!(found_f && found_i))
        {
          std::cerr << "Error performing qsq sort, need to include rest matrix elements" << std::endl;
      
          if(!!!found_f)
            std::cerr << "Couldn't find the final state" << std::endl;

          if(!!!found_i)
            std::cerr << "Couldn't find the initial state" << std::endl;


          dump_xml(); 
          exit(1); // abort and write out bad xml files
        }

        const double xi = unsorted.begin()->m_ini.xi; 
        const double L_s = unsorted.begin()->m_ini.L_s; 
        const double factor = 1./xi * 2. * acos(-1.) / L_s; 

        std::map<double,std::vector<CartesianMatrixElement> > m_map; 
        std::map<double,std::vector<CartesianMatrixElement> >::iterator mapit; 
        for(it = unsorted.begin(); it != unsorted.end(); ++it)
        {
          // skip the baddies
          if(!!!it->have_active_data)
            continue; 

          // ignore any type of level splittings across different irreps 
          // also ignore any statistical fluctuations associated w/ different mom directions
          double q2 = Mink_qsq(*it, E_f, E_i, factor);           
          mapit = m_map.find(q2); 

          if(mapit == m_map.end()) 
          {
            std::vector<CartesianMatrixElement> dum(1,*it); 
            m_map.insert(std::pair<double,std::vector<CartesianMatrixElement> >(q2,dum)); 
          }
          else
            mapit->second.push_back(*it); 

        }

        sorted_by_Q2.clear(); 
        for(mapit = m_map.begin(); mapit != m_map.end(); ++mapit)
          sorted_by_Q2.push_back(mapit->second); 
      }



      // this did not actially work.. 

#if 0
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
        for(map_it = my_map.begin(); map_it != my_map.end(); ++map_it)
          sorted_by_Q2.push_back(map_it->second); 

      }

#endif 


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


      ENSEM::EnsemReal computeQ2(const std::vector<CartesianMatrixElement> &elems)
      {
        std::vector<CartesianMatrixElement>::const_iterator it;
        for(it = elems.begin(); it != elems.end(); ++it)
          if(it->have_active_data)
            break;

        POW2_ASSERT(it != elems.end()); 

        ENSEM::EnsemReal E_f , E_i, Q;
        ADATXML::Array<int> p_f,p_i,q;
        double xi = it->m_ini.xi;
        double L_s = it->m_ini.L_s;
        double mom_fac =  1./xi * 2.*acos(-1.)/L_s;
        int qq; 
        E_f = it->E_f;
        E_i = it->E_i;
        p_f = it->p_f;
        p_i = it->p_i; 

        Q = E_i - E_f; 
        q = p_i;
        q[0] -= p_f[0];
        q[1] -= p_f[1];
        q[2] -= p_f[2];
        qq = 0; 
        qq += q[0]*q[0];
        qq += q[1]*q[1];
        qq += q[2]*q[2];

        return (-Q*Q + SEMBLE::toScalar(mom_fac*mom_fac*double(qq)));
      }


      void write_projected_matrix_elements(void)
      {
        std::vector<std::vector<CartesianMatrixElement> >::const_iterator it; 

        for(it = sorted_by_Q2.begin(); it != sorted_by_Q2.end(); ++it)
          write_proj(*it); 
      }


      void write_proj(const std::vector<CartesianMatrixElement> &elems)
      {
        double Q2 = SEMBLE::toScalar(ENSEM::mean(computeQ2(elems))); 
        std::vector<CartesianMatrixElement>::const_iterator it;
        for(it = elems.begin(); it != elems.end(); ++it)
          write_proj(*it, Q2);
      }

      void write_proj(const CartesianMatrixElement &elem, const double q2)
      {
        if(!!!elem.have_active_data)
          return;

        std::string path = SEMBLE::SEMBLEIO::getPath();
        std::stringstream ss;
        ss << path << "cont_projected_matrix_elements";
        SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
        ss << "/Q2_" << q2;
        SEMBLE::SEMBLEIO::makeDirectoryPath(ss.str());
        write_proj(elem.t,elem.m_elem,ss.str(),"t");
        write_proj(elem.x,elem.m_elem,ss.str(),"x");
        write_proj(elem.y,elem.m_elem,ss.str(),"y");
        write_proj(elem.z,elem.m_elem,ss.str(),"z");
      }

      std::string elem_str(const simpleWorld::ContinuumStatePrimitive &s)
      {
        std::stringstream ss;
        ss << "J" << s.J << "_H" << s.H << "_p" << s.mom[0] << s.mom[1] << s.mom[2]; 
        return ss.str(); 
      }

      void write_proj(const accumulator &a, const simpleWorld::ContinuumMatElem &elem,
          const std::string &pth, const std::string &comp)
      {
        if(!!!a.active)
          return;
        if(!!!a.success)
          return;

        std::stringstream fname; 
        fname << pth << "/SINK_V_SOURCE_" << elem_str(elem.sink.state) 
          << "__V_" << comp << "__" << elem_str(elem.source.state) ;

        std::string real, imag, jack;
        real = fname.str() + std::string("__real.ax");
        imag = fname.str() + std::string("__imag.ax");
        jack = fname.str() + std::string("__corr.jack"); 

        AxisPlot preal, pimag; 

        preal.addEnsemData(ENSEM::real(a.corr),"//sq",1);
        pimag.addEnsemData(ENSEM::imag(a.corr),"//sq",1);

        preal.sendToFile(real);
        pimag.sendToFile(imag);
        ENSEM::write(jack,a.corr); 

      }





      // data 
      ThreePointCorrIni_t m_ini; 
      std::vector<std::vector<CartesianMatrixElement> > sorted_by_Q2; 
      std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> m_bad_norms; 

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

        ENSEM::EnsemReal Q2 = work_manager.computeQ2(*it_q2); 
        tmp->setQ2(Q2);
        tmp->zeroFilter(); 
        ret.push_back(tmp);
      }

      return ret; 
    }


}




