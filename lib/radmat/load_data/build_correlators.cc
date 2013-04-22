/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : build_correlators.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Mon Apr 22 09:09:51 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "build_correlators.h"
#include "invert_subduction.h"
#include "simple_world.h"
#include "generate_redstar_xml.h"  
#include "radmat_database_interface.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/perThreadStorage.h"
#include "radmat/llsq/llsq_gen_system.h"
#include "semble/semble_semble.h"
#include "radmat_overlap_key_val_db.h"
#include "hadron/ensem_filenames.h"
#include "io/adat_xmlio.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"
#include "jackFitter/plot.h"
#include "ensem/ensem.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <omp.h>



// #define DEBUG_CORRELATOR_NORMALIZATION // do loads of printing at the normalization stage
// #define SERIOUSLY_DEBUG_CORRELATOR_NORMALIZATION  // turn on annoying printing

// #define BUILD_CORRS_USE_OMP_PARALLEL

// #define  DEBUG_AT_MAKE_MOM_INV_TAGS // are we making the tags right?


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
      << " isProjected = " << o.isProjected << " maSource = " << o.maSource
      << " maSink = " << o.maSink; 
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
    doXMLRead(ptop,"maSource",prop.maSource,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"maSink",prop.maSink,__PRETTY_FUNCTION__); 
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
    write(xml,"maSource",prop.maSource);
    write(xml,"maSink",prop.maSink); 
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
        {
#ifdef DEBUG_CORRELATOR_NORMALIZATION
          std::cout << __func__ << ": unsuccessful " << std::endl;
#endif
          return;
        }

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
        momentum_factor = mom_factor(ini.xi,ini.L_s);      
        p_f = it->m_obj.redstar_xml.npoint[1].irrep.mom; 
        p_i = it->m_obj.redstar_xml.npoint[3].irrep.mom;

        RadmatMassOverlapData_t source_tmp = db.fetch(it->m_obj.source_normalization); 

        E_f = source_tmp.E()*SEMBLE::toScalar(0.);
        E_i = E_f;
        int ct = 0;

        for(it = data.begin(); it != data.end(); ++it)
        {
          ENSEM::EnsemVectorComplex corr_tmp = db.fetch(it->m_obj.redstar_xml);
          RadmatMassOverlapData_t source = db.fetch(it->m_obj.source_normalization); 
          RadmatMassOverlapData_t sink = db.fetch(it->m_obj.sink_normalization); 

#ifdef DEBUG_CORRELATOR_NORMALIZATION
          std::string outstem = Hadron::ensemFileName(it->m_obj.redstar_xml); 
          std::string pth = SEMBLE::SEMBLEIO::getPath();
          std::stringstream path;
          path << pth << "debug_normalization";
          SEMBLE::SEMBLEIO::makeDirectoryPath(path.str());
          path << "/" << outstem;
          ENSEM::write(path.str() + std::string("_corr_pre") , corr_tmp); 
          ENSEM::EnsemVectorComplex norm = corr_tmp * SEMBLE::toScalar(0.); 
#endif

          // the hadron key uses 1 based arrays

          // NB: assumption that npt is organized like <sink, ins , source>

          const int t_source(it->m_obj.redstar_xml.npoint[3].t_slice);
          const int t_sink(it->m_obj.redstar_xml.npoint[1].t_slice); 


          POW2_ASSERT(t_source < t_sink); 

#ifdef  SERIOUSLY_DEBUG_CORRELATOR_NORMALIZATION 
          std::cout << __func__ << std::endl;
          std::cout << "Z_sink = " << SEMBLE::toScalar(ENSEM::mean(sink.Z())) << std::endl;
          std::cout << "Z_source = " << SEMBLE::toScalar(ENSEM::mean(source.Z())) << std::endl;
          std::cout << "E_sink = " << SEMBLE::toScalar(ENSEM::mean(sink.E())) << std::endl;
          std::cout << "E_source = " << SEMBLE::toScalar(ENSEM::mean(source.E())) << std::endl;
          std::cout << "tsink = " << t_sink << " tsource =   " << t_source << std::endl;

#endif 

          

          // NB: the indexing here assumes [tsource,tsink] ie: inclusive range
          for(int t_ins = t_source; t_ins <= t_sink; ++t_ins)
          {

            ENSEM::EnsemReal prop = propagation_factor(sink.E(),sink.Z(),t_sink,t_ins,
                source.E(),source.Z(),t_source);

            ENSEM::pokeObs(corr_tmp,ENSEM::peekObs(corr_tmp,t_ins)/prop,t_ins);

#ifdef SERIOUSLY_DEBUG_CORRELATOR_NORMALIZATION
            if(t_ins == 0) 
              std::cout << "t_ins = " << t_ins << " norm = " << SEMBLE::toScalar(ENSEM::mean(prop)) << std::endl;
#endif



#ifdef DEBUG_CORRELATOR_NORMALIZATION

            ENSEM::pokeObs(norm,prop,t_ins); 
#endif
          } // end loop over t_ins

#ifdef DEBUG_CORRELATOR_NORMALIZATION
          ENSEM::write(path.str() + std::string("_corr_post"), corr_tmp);
          ENSEM::write(path.str() + std::string("_norm") , norm);
  
          std::cout << __func__ << std::endl; 
          std::cout << SEMBLE::toScalar(it->m_coeff) << " X " << Hadron::ensemFileName(it->m_obj.redstar_xml) << std::endl;
#endif      


          E_f = E_f + sink.E();
          E_i = E_i + source.E();

          corr = corr + it->m_coeff*corr_tmp; 

          ++ct;

        }  // end loop over data (iterator loop)

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


      CartesianMatrixElement(const ThreePointCorrIni_t &ini, 
          const simpleWorld::ContinuumMatElem &elem)
        : have_active_data(false) , success(true) , m_ini(ini) , m_elem(elem) 
      {

        am_i = ini.threePointCorrXMLIni.maSource;
        am_f = ini.threePointCorrXMLIni.maSink;
        // hardwire some assumptions about projected 
        // operator naming here.. if something changes we can just build some
        // sort of hash to retreive the correct name but the current
        // convetions should be good for a bit
        if(m_ini.threePointCorrXMLIni.isProjected)
        {
          m_elem.source.state.name = append_momentum(m_elem.source.state.name,
              m_elem.source.state.mom);
          m_elem.sink.state.name = append_momentum(m_elem.sink.state.name,
              m_elem.sink.state.mom);
        }

        if(m_ini.threePointCorrXMLIni.isDiagonal)
          if(am_i != am_f)
          {
            std::cerr << __func__ << ": WARNING: you chose a diagonal matrix element"
              <<" but gave different rest masses for source and sink. " << std::endl;
          }

        redstarCartMatElem celem(m_elem,ini.threePointCorrXMLIni.source_id,
            ini.threePointCorrXMLIni.sink_id);
        t.eat(celem.get_component(0),m_ini);
        x.eat(celem.get_component(1),m_ini);
        y.eat(celem.get_component(2),m_ini);
        z.eat(celem.get_component(3),m_ini);

        check(); 

        std::stringstream ss; 
        ss << m_ini.matElemID << "_" << m_elem.sink.state.H << "_" << m_elem.source.state.H;
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
        else
          have_active_data = false;

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
        ret.mom_fac = mom_factor(m_ini.xi , m_ini.L_s);     
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


      std::string write_me(void) const
      {
        std::stringstream ss; 
        ss << m_mat_elem_id << ".p_f"  << p_f[0] << p_f[1] << p_f[2] << ".p_i"  
          << p_i[0] << p_i[1] << p_i[2] << ".am_f" << std::setw(3) << am_f << ".am_i" 
          << std::setw(3) << am_i; 
        return ss.str();  
      }


      bool have_active_data;
      bool success;
      double am_f, am_i; 
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



    // this is the main work horse..
    std::vector<CartesianMatrixElement> 
      getCartesianMatrixElements(const std::vector<simpleWorld::ContinuumMatElem> & elems, const ThreePointCorrIni_t &ini)
      {
        std::vector<CartesianMatrixElement> ret;
        std::vector<simpleWorld::ContinuumMatElem>::const_iterator it; 
        for(it = elems.begin(); it != elems.end(); ++it)
          ret.push_back(CartesianMatrixElement(ini,*it)); 
        return ret; 
      }


    // this is a threaded main workhorse -- on my mac(4core, 8threads) I saw a factor of ~2 speed up running parallel 
    std::vector<CartesianMatrixElement>
      getCartesianMatrixElementsParallel(const std::vector<simpleWorld::ContinuumMatElem> &elems, const ThreePointCorrIni_t &ini)
      {
        /*
           std::cout << __func__ << ": debugging is on, elems contains" << std::endl;
           std::vector<simpleWorld::ContinuumMatElem>::const_iterator it; 
           for(it = elems.begin(); it != elems.end(); ++it)
           std::cout << *it << std::endl;
         */

        std::vector<CartesianMatrixElement> collect;
        registerSubductionTables(); // do this mono to avoid race 
        int index;
        int sz = elems.size();
        collect.resize(sz);

#pragma omp parallel for shared(index,sz)

        for(index =0; index < sz; ++index)
          collect[index] = CartesianMatrixElement(ini,elems[index]);

        return collect; 
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


    double Mink_qsq(const CartesianMatrixElement & elem, 
        const double E_f , const double E_i, const double factor)
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

#ifdef  BUILD_CORRS_USE_OMP_PARALLEL 

        unsorted = getCartesianMatrixElementsParallel(
            simpleWorld::getContinuumMatElemFromXML(m_ini.threePointCorrXMLIni.continuumMatElemXML),
            m_ini);

#else

        unsorted = getCartesianMatrixElements(
            simpleWorld::getContinuumMatElemFromXML(m_ini.threePointCorrXMLIni.continuumMatElemXML),
            m_ini);

#endif

        sort(unsorted,m_ini.threePointCorrXMLIni.isDiagonal);

        dump_unsorted(unsorted);
        dump_sorted();

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


        //     std::cout << __func__ << ": unsorted.size() " << unsorted.size() << std::endl;


        double E_f, E_i; // rest energies..       
        std::vector<CartesianMatrixElement>::const_iterator it;

        E_f = unsorted.begin()->am_f; 
        E_i = unsorted.begin()->am_i;


        for(it = unsorted.begin(); it != unsorted.end(); ++it)
        {

          if(!!!it->success)
            push_bad_list(*it);

        }


        const double xi = unsorted.begin()->m_ini.xi; 
        const double L_s = unsorted.begin()->m_ini.L_s; 
        const double factor = mom_factor(xi , L_s); 

        std::map<double,std::vector<CartesianMatrixElement> > m_map; 
        std::map<double,std::vector<CartesianMatrixElement> >::iterator mapit; 
        for(it = unsorted.begin(); it != unsorted.end(); ++it)
        {
          // skip the baddies
          if(!!!it->have_active_data)
          {
            std::cout << __func__ << ": skipping a baddie" << std::endl;
            continue; 
          }
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


        //      std::cout << __func__ << ": sorted_by_Q2.size() = " << sorted_by_Q2.size() << std::endl;
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
        double mom_fac =  mom_factor(xi , L_s);
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




      void dump_unsorted(const std::vector<CartesianMatrixElement> &elems)
      {
        std::ofstream out("unsorted_cart.txt");
        std::vector<CartesianMatrixElement>::const_iterator it;
        for(it = elems.begin(); it != elems.end(); ++it)
          out << it->write_me() << "\n";
        out.close(); 
      }


      void dump_sorted(void)
      {
        std::vector<std::vector<CartesianMatrixElement> >::const_iterator it; 
        std::vector<CartesianMatrixElement>::const_iterator vit;
        std::ofstream out("sorted_cart.txt");

        for(it = sorted_by_Q2.begin(); it != sorted_by_Q2.end(); ++it)
        {
          out << "MARK////////////////////////////" << std::endl;
          for(vit = it->begin(); vit != it->end(); ++vit)
            out << vit->write_me() << std::endl; 
        }

        out.close(); 

      }



      // data 
      ThreePointCorrIni_t m_ini; 
      std::vector<std::vector<CartesianMatrixElement> > sorted_by_Q2; 
      std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
      std::vector<RadmatExtendedKeyHadronNPartIrrep_t> m_bad_norms; 

    }; // manageWork


    LatticeMultiDataTag get_lattice_multi_data_tag(const CartesianMatrixElement &c, const int jmu)
    {
      LatticeMultiDataTag out;

      std::stringstream ss;
      ss << simpleWorld::stateFileName(c.m_elem.sink) << ".V_" << jmu << "." << simpleWorld::stateFileName(c.m_elem.source); 
      out.file_id = ss.str(); 
      out.jmu = jmu;
      out.mat_elem_id = c.m_mat_elem_id;
      out.p_f = c.p_f; 
      out.p_i = c.p_i;
      out.E_f = c.E_f;
      out.E_i = c.E_i;
      out.mom_fac = mom_factor(c.m_ini.xi, c.m_ini.L_s);




#ifdef DEBUG_AT_MAKE_MOM_INV_TAGS 

      std::cout << __func__ << ": debuggin on" << std::endl;
      out.print_me();
      std::cout << "mom_string = " << out.mom_string() << std::endl;
      std::cout << "E_string = " << out.E_string() << std::endl;
#endif



      return out;
    }


    ADAT::Handle<LLSQLatticeMultiData> get_lattice_multi_data(const std::vector<CartesianMatrixElement> &in)
    {
      ADAT::Handle<LLSQLatticeMultiData> out(new LLSQLatticeMultiData);
      POW2_ASSERT(&*out); // alloc check 
      std::vector<CartesianMatrixElement>::const_iterator it; 


      for(it = in.begin(); it != in.end(); ++it)
      {
        if(!!!it->have_active_data)
          continue; 

        LatticeMultiDataTag t,x,y,z; 
        bool tt,xx,yy,zz;     

        tt = it->t.success && it->t.active; 
        xx = it->x.success && it->x.active;
        yy = it->y.success && it->y.active;
        zz = it->z.success && it->z.active; 


        if(tt)
        {
          //          SEMBLE::SembleVector<std::complex<double> > foo;
          //          foo = it->t.corr; 
          //          std::cout << "t = " << foo.mean() << std::endl;
          out->append_row_ensem(it->t.corr,get_lattice_multi_data_tag(*it,0));
        } 
        if(xx)
        {
          //          SEMBLE::SembleVector<std::complex<double> > foo;
          //          foo = it->x.corr; 
          //          std::cout << "x = " << foo.mean() << std::endl;
          out->append_row_ensem(it->x.corr,get_lattice_multi_data_tag(*it,1));
        }
        if(yy)
        {
          //          SEMBLE::SembleVector<std::complex<double> > foo;
          //          foo = it->y.corr; 
          //          std::cout << "y = " << foo.mean() << std::endl;
          out->append_row_ensem(it->y.corr,get_lattice_multi_data_tag(*it,2));
        }
        if(zz)
        {          
          //          SEMBLE::SembleVector<std::complex<double> > foo;
          //          foo = it->z.corr; 
          //          std::cout << "z = " << foo.mean() << std::endl;
          out->append_row_ensem(it->z.corr,get_lattice_multi_data_tag(*it,3)); 
        }

      }


      return out; 

    }









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

          /*
             std::vector<LLSQDataPoint>::const_iterator print_it; 
             std::cout << __func__ << ": t_ins = " << t_ins << std::endl;
             for(print_it = points_at_this_t_ins.begin(); print_it != points_at_this_t_ins.end(); ++print_it)
             std::cout << print_it->toString() << std::endl;
           */

          tmp->insert(LLSQDataPointQ2Pack::value_type(t_ins,points_at_this_t_ins)); 
        }

        ENSEM::EnsemReal Q2 = work_manager.computeQ2(*it_q2); 
        tmp->setQ2(Q2);
        tmp->zeroFilter();

        // only send back guys with active data into the llsq
        if(!!!tmp->haveData())
        {
          std::cout << __func__ << ": removing Q2 = " << SEMBLE::toScalar(ENSEM::mean(Q2)) 
            << " from llsq system (no active data) " << std::endl;
          continue; 
        }
        else
          ret.push_back(tmp);
      }


      return ret; 
    }


  std::vector<ADAT::Handle<LLSQLatticeMultiData> > BuildCorrelators::build_multi_correlators(void)
  {
    POW2_ASSERT(have_ini);
    manageWork work_manager(m_ini);
    std::vector<std::vector<CartesianMatrixElement> > sorted_by_q2(work_manager.get_sorted_data());
    std::vector<std::vector<CartesianMatrixElement> >::const_iterator it_q2;
    std::vector<ADAT::Handle<LLSQLatticeMultiData> >  ret; 

    // this is looping over all of the Q2 values
    for(it_q2 = sorted_by_q2.begin(); it_q2 != sorted_by_q2.end(); ++it_q2)
      ret.push_back(get_lattice_multi_data(*it_q2)); 

    return ret;
  }


}


#undef BUILD_CORRS_USE_OMP_PARALLEL
#undef DEBUG_CORRELATOR_NORMALIZATION 
#undef SERIOUSLY_DEBUG_CORRELATOR_NORMALIZATION  
#undef DEBUG_AT_MAKE_MOM_INV_TAGS 



