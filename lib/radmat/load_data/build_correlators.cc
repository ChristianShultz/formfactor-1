/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : build_correlators.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Fri 01 Nov 2013 04:07:49 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "build_correlators.h"
#include "build_correlators_bad_data_repository.h"
#include "lattice_multi_data_tag.h"
#include "lattice_multi_data_object.h"
#include "invert_subduction.h"
#include "g_parity_world.h"
#include "g_parity_world_generate_redstar_xml.h"  
#include "radmat_database_interface.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/perThreadStorage.h"
#include "radmat/utils/mink_qsq.h"
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
#include <math.h>
#include <omp.h>


 #define BUILD_CORRELATORS_PARALLEL

// #define DEBUG_BUILD_LLSQ_DATA

namespace radmat
{




  // BUILD CORRELATORS
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////




  typedef radmatAllConfDatabaseInterface< Hadron::KeyHadronNPartNPtCorr_t,
          ENSEM::EnsemVectorComplex,
          RadmatExtendedKeyHadronNPartIrrep_t,
          RadmatMassOverlapData_t> DatabaseInterface_t;


  // utility
  namespace
  {


    LatticeMultiDataTag get_lattice_tag(const int jmu, 
        const ThreePointCorrIni_t &ini, 
        const gParityWorld::GParityHelicityMatrixElement &e)
    {
      LatticeMultiDataTag ret; 

      std::stringstream ss;
      ss << gParityWorld::stateFileName(e.sink)  << ".V_" << jmu 
        << "." << gParityWorld::stateFileName(e.source);

      ret.file_id = ss.str();
      ret.jmu = jmu; 

      std::stringstream mat_elem_id;
      mat_elem_id << ini.matElemID << "_" << e.sink.state.H 
        << "_" << e.source.state.H; 

      ret.mat_elem_id = mat_elem_id.str(); 
      ret.p_f = e.sink.state.mom; 
      ret.p_i = e.source.state.mom;

      // a sanity check, we will always assum a size of 3 
      if(ret.p_f.size() != 3 || ret.p_i.size() != 3)
      {
        std::cerr << __func__ << ": sanity check failed here" << std::endl;
        exit(1); 
      }

      ret.mom_fac = mom_factor(ini.xi,ini.L_s); 
      ret.set_qsq_label(Mink_qsq(ret.p_f,ini.threePointCorrXMLIni.maSink,
            ret.p_i, ini.threePointCorrXMLIni.maSource,ret.mom_fac)); 

      return ret; 
    }


    struct cartesianMatrixElementXML
    {        
      typedef generateCircularRedstarXML::data_t data_t;
      typedef mergeSubducedOperators3pt::NPointKey NPointKey;
      typedef mergeSubducedOperators3pt::objNPointKey objNPointKey; 
      typedef mergeSubducedOperators3pt::listNPointKey listNPointKey;

      cartesianMatrixElementXML(const int lorentz_index, 
          const data_t &d, 
          const gParityWorld::GParityHelicityMatrixElement &e,
          const ThreePointCorrIni_t &ini)
        :have_active_data(d.first), my_npoint(d.second) , my_helicity_elem(e) 
      {
        my_tag = get_lattice_tag(lorentz_index,ini,e);
      }

      double qsq_tag(void) const {return my_tag.get_qsq_label();}

      bool have_active_data; 
      LatticeMultiDataTag my_tag;
      listNPointKey my_npoint; 
      gParityWorld::GParityHelicityMatrixElement my_helicity_elem; 
    };

    std::ostream & operator<<(std::ostream & o , const cartesianMatrixElementXML &e)
    {
      o << e.my_helicity_elem;
      return o;
    }


    // brute force sort on qsq
    std::map<double,std::vector<cartesianMatrixElementXML> > 
      sortByQ2(const std::vector<cartesianMatrixElementXML> &unsorted)
      {
        std::map<double,std::vector<cartesianMatrixElementXML> > ret; 
        std::vector<cartesianMatrixElementXML>::const_iterator it; 

        // loop everything and toss it into a vector if it matches else make a new vector
        for(it = unsorted.begin(); it != unsorted.end(); ++it)
        {

          if(ret.find(it->qsq_tag()) != ret.end())
            ret.find(it->qsq_tag())->second.push_back(*it);
          else
            ret.insert(std::map<double,std::vector<cartesianMatrixElementXML> >::value_type(
                  it->qsq_tag(),std::vector<cartesianMatrixElementXML>(1,*it))); 
        }

        std::map<double,std::vector<cartesianMatrixElementXML> >::const_iterator it1; 
        int ct(0);
        for(it1 = ret.begin(); it1 != ret.end(); ++it1)
          ct += it1->second.size();

        return ret;
      }


/*
 *  A LOCAL bad data repository that exists only here
 *
 */

    BuildCorrsLocalBadDataRepo_t bad_data_repository;

    // NB hardwire in here
    struct dbKeyObject
    {    
      typedef mergeSubducedOperators3pt::NPointKey nptkey;
      typedef RadmatExtendedKeyHadronNPartIrrep_t normkey;

      dbKeyObject(const nptkey &k , const std::string &id_sink, const std::string &id_source)
        : npt(k)
      {
        sink = normkey(id_sink,k.npoint[1].irrep);
        source = normkey(id_source,k.npoint[3].irrep); 
      }

      nptkey npt;
      normkey sink,source;  
    };  

    typedef ObjExpr_t<ENSEM::Complex,dbKeyObject> coeffKey;
    typedef ListObjExpr_t<ENSEM::Complex,dbKeyObject> listCoeffKey; 

    // get the cartesian matrix element into a simple form 
    listCoeffKey genListCoeffKey(const cartesianMatrixElementXML &e, 
        const std::string &sink, 
        const std::string &source)
    {
      listCoeffKey ret; 
      cartesianMatrixElementXML::listNPointKey::const_iterator it;
      for(it = e.my_npoint.begin(); it != e.my_npoint.end(); ++it)
        ret = ret + coeffKey(it->m_coeff, dbKeyObject(it->m_obj,sink,source));
      return ret; 
    }





    struct build_correlator_return_value
    {
      LatticeMultiDataTag tag;
      ENSEM::EnsemVectorComplex data;
      bool success; 
    };




    build_correlator_return_value
      build_a_correlator( std::vector<Hadron::KeyHadronNPartNPtCorr_t >  & bad_corrs,
          std::vector<RadmatExtendedKeyHadronNPartIrrep_t> & bad_norms,
          const cartesianMatrixElementXML &d,
          const ThreePointCorrIni_t &ini,
          const DatabaseInterface_t &db)
      {
        build_correlator_return_value ret; 
        ret.tag = d.my_tag;

        std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
        std::vector<RadmatExtendedKeyHadronNPartIrrep_t> m_bad_norms; 
        listCoeffKey data = 
          genListCoeffKey(d,ini.threePointCorrXMLIni.sink_id,ini.threePointCorrXMLIni.source_id); 
        listCoeffKey::const_iterator it;
        bool success = true; 

#ifdef DEBUG_BUILD_LLSQ_DATA
        int counter = 1; 
#endif 

        for(it = data.begin(); it != data.end(); ++it)
        {
          if(!!!db.exists(it->m_obj.npt))
            m_bad_corrs.push_back(it->m_obj.npt);
          if(!!!db.exists(it->m_obj.source))
            m_bad_norms.push_back(it->m_obj.source);
          if(!!!db.exists(it->m_obj.sink))
            m_bad_norms.push_back(it->m_obj.sink);

#ifdef DEBUG_BUILD_LLSQ_DATA 
          std::cout << __func__ << " check # " << counter << " of " 
            << data.size() << " checked for " << it->m_obj.npt << std::endl; 
          ++counter;
#endif
        }

        // break early if not successful 
        if(!!!m_bad_corrs.empty())
          success = false; 

        if(!!!m_bad_norms.empty())
          success = false; 

        // this shouldn't happen but would 
        // probably cause some nasty if it did
        if(data.begin() == data.end())
        {
          std::cerr << __func__ << ": somehow you pumped in an empty guy" 
            << " I will try to continue" << std::endl;
          success = false; 
        }

        ret.success = success;

        // attempt a graceful exit
        if(!!!success)
        {
          // set up a zero valued three point correlator
          ENSEM::EnsemVectorComplex corr; 
          corr.resize(1);
          corr.resizeObs(1);
          corr = SEMBLE::toScalar(0.); 


          bad_corrs.insert(bad_corrs.end(),m_bad_corrs.begin(),m_bad_corrs.end());
          bad_norms.insert(bad_norms.end(),m_bad_norms.begin(),m_bad_norms.end());
          ret.data = corr; 
          return ret;
        }

        // set up a zero valued three point correlator
        ENSEM::EnsemVectorComplex corr; 
        it = data.begin();
        corr = db.fetch(it->m_obj.npt);
        corr = SEMBLE::toScalar(0.); 

        // initialize some variables        
        ThreePtPropagationFactor<double> propagation_factor;
        ADATXML::Array<int> p_f = it->m_obj.npt.npoint[1].irrep.mom; 
        ADATXML::Array<int> p_i = it->m_obj.npt.npoint[3].irrep.mom;

        RadmatMassOverlapData_t source_tmp = db.fetch(it->m_obj.source); 

        ENSEM::EnsemReal E_f = source_tmp.E()*SEMBLE::toScalar(0.);
        ENSEM::EnsemReal E_i = E_f;
        int ct = 0;



        for(it = data.begin(); it != data.end(); ++it)
        {
#ifdef DEBUG_BUILD_LLSQ_DATA 
          std::cout << "\n\n" << __func__ << ": BUILDING " << it->m_obj.npt << std::endl; 
#endif
          ENSEM::EnsemVectorComplex corr_tmp = db.fetch(it->m_obj.npt);
          RadmatMassOverlapData_t source = db.fetch(it->m_obj.source); 
          RadmatMassOverlapData_t sink = db.fetch(it->m_obj.sink); 

          std::string outstem = Hadron::ensemFileName(it->m_obj.npt); 
          std::string pth = SEMBLE::SEMBLEIO::getPath();
          std::stringstream path;

          path << pth << "Q2_" << ret.tag.get_qsq_label(); 
          SEMBLE::SEMBLEIO::makeDirectoryPath(path.str()); 

          path << "/correlator_normalization";
          SEMBLE::SEMBLEIO::makeDirectoryPath(path.str());

          path << "/" << outstem;
          ENSEM::write(path.str() + std::string("_corr_pre") , corr_tmp); 
          ENSEM::EnsemVectorComplex norm = corr_tmp * SEMBLE::toScalar(0.); 

          // the hadron key uses 1 based arrays

          // NB: assumption that npt is organized like <sink, ins , source>
          const int t_source(it->m_obj.npt.npoint[3].t_slice);
          const int t_sink(it->m_obj.npt.npoint[1].t_slice); 

          // sanity
          POW2_ASSERT(t_source < t_sink); 


          // NB: the indexing here assumes [tsource,tsink] ie: inclusive range
          for(int t_ins = t_source; t_ins <= t_sink; ++t_ins)
          {

            ENSEM::EnsemReal prop = propagation_factor(sink.E(),sink.Z(),t_sink,t_ins,
                source.E(),source.Z(),t_source);

            ENSEM::pokeObs(corr_tmp,ENSEM::peekObs(corr_tmp,t_ins)/prop,t_ins);

            ENSEM::pokeObs(norm,prop,t_ins); 
          } // end loop over t_ins

          ENSEM::write(path.str() + std::string("_corr_post"), corr_tmp);
          ENSEM::write(path.str() + std::string("_norm") , norm);


          E_f = E_f + sink.E();
          E_i = E_i + source.E();

          corr = corr + it->m_coeff*corr_tmp; 

          ++ct;

        }  // end loop over data (iterator loop)

        // --- DO NORMALIZATION/RENORMALIZATION HERE
        // these are real but the * operator is only overloaded
        // for complex types so shove a zero down its throat
        // to get the desired result
        double Z_t_d , Z_s_d ; 

        // the wonderful flavor wick contraction factor
        Z_t_d = ini.threePointCorrXMLIni.renormalization.RGE_prop; 
        Z_s_d = Z_t_d; 

        // the thing you would quote in a paper
        Z_t_d *= ini.threePointCorrXMLIni.renormalization.Z_t; 
        Z_s_d *= ini.threePointCorrXMLIni.renormalization.Z_s; 
        ENSEM::Complex Z_t, Z_s; 

        Z_t = SEMBLE::toScalar(std::complex<double>(Z_t_d,0.)); 
        Z_s = SEMBLE::toScalar(std::complex<double>(Z_s_d,0.));  

        // fill in the rest of the relevant data
        E_f = E_f/SEMBLE::toScalar(double(ct));
        E_i = E_i/SEMBLE::toScalar(double(ct));

        if ( ret.tag.mu() == 0 ) 
          ret.data = corr / Z_t;
        else
          ret.data = corr / Z_s;  

        ret.tag.E_f = E_f;
        ret.tag.E_i = E_i;


#ifdef DEBUG_BUILD_LLSQ_DATA 
        std::cout << __func__ << ": " << ret.tag.splash_tag() << std::endl;
#endif

        return ret;
      };




    std::pair<bool, ADAT::Handle<LLSQLatticeMultiData> >
      build_llsq_data(const std::vector<cartesianMatrixElementXML> &e,
          const ThreePointCorrIni_t &ini)
      {
        ADAT::Handle<LLSQLatticeMultiData> ret(new LLSQLatticeMultiData()); 
        std::vector<cartesianMatrixElementXML>::const_iterator it;
        std::string idsource,idsink;
        idsink = ini.threePointCorrXMLIni.sink_id;
        idsource = ini.threePointCorrXMLIni.source_id; 

        DatabaseInterface_t db(ini.radmatDBProp); 
        std::vector<Hadron::KeyHadronNPartNPtCorr_t > m_bad_corrs;
        std::vector<RadmatExtendedKeyHadronNPartIrrep_t> m_bad_norms; 

        bool any_data = false; 

        for(it = e.begin(); it != e.end(); ++it)
        {
          if(!!!it->have_active_data) 
            continue; 
          build_correlator_return_value tmp;
          tmp = build_a_correlator(m_bad_corrs,
              m_bad_norms,
              *it,
              ini,
              db
              );

          if(tmp.success)
            ret->append_row_ensem(tmp.data,tmp.tag);

          if(tmp.success)
            any_data = true; 
        }

        any_data &= !!!ret->empty(); 

        bad_data_repository.insert(omp_get_thread_num(),m_bad_corrs); 
        bad_data_repository.insert(omp_get_thread_num(),m_bad_norms); 

        return std::pair<bool,ADAT::Handle<LLSQLatticeMultiData> >(any_data,ret); 
      }


    std::vector<ADAT::Handle<LLSQLatticeMultiData> >
      build_llsq_data_parallel(const std::map<double, std::vector<cartesianMatrixElementXML> > &m,
          const ThreePointCorrIni_t &ini)
      {
        std::vector<bool> use_data(m.size(),false); 
        std::vector<ADAT::Handle<LLSQLatticeMultiData> > tmp,ret; 
        tmp.resize(m.size()); 
        std::vector<std::pair<double,std::vector<cartesianMatrixElementXML> > > elems;
        std::map<double,std::vector<cartesianMatrixElementXML> >::const_iterator it; 

        // populate a vector for parallel
        for(it = m.begin(); it != m.end(); ++it)
          elems.push_back(
              std::pair<double,std::vector<cartesianMatrixElementXML> >(it->first,it->second)); 

        int elem;
#pragma omp parallel for shared(elem)  schedule(dynamic,1)

        for(elem = 0; elem < elems.size(); ++elem)
        {
          //         std::cout << "Working on Q2 = " << elems[elem].first << std::endl; 
          std::pair<bool,ADAT::Handle<LLSQLatticeMultiData> > val; 
          val = build_llsq_data(elems[elem].second,ini);

          use_data[elem] = val.first;
          tmp[elem] = val.second;
        }

#pragma omp barrier 

        for(elem = 0; elem < elems.size(); ++elem)
        {
          if(use_data[elem])
          {
            ret.push_back(tmp[elem]); 
          }
          else
          {
            std::cout << "removing Q2 = " << elems[elem].first 
              << " from llsq (no data)" << std::endl;
          }
        }
        return ret; 
      }






  } // anonymous 




  std::vector<ADAT::Handle<LLSQLatticeMultiData> >
    BuildCorrelators::build_multi_correlators(void)
    {
      POW2_ASSERT(have_ini); 

      // need thes
      registerSubductionTables(); 

      // transform the xml into the internal helicity representation
      std::vector<gParityWorld::GParityHelicityMatrixElement> internal
        = gParityWorld::getGParityHelicityMatrixElementFromXML(
            m_ini.threePointCorrXMLIni.continuumMatElemXML); 

      std::vector<gParityWorld::GParityHelicityMatrixElement>::const_iterator internal_it;
      std::vector<cartesianMatrixElementXML> unsorted_elements; 
      std::vector<cartesianMatrixElementXML>::const_iterator unsorted_iterator; 
      std::map<double, std::vector<cartesianMatrixElementXML> > sorted_elements;
      std::map<double, std::vector<cartesianMatrixElementXML> >::const_iterator sorted_it;  
      std::vector<ADAT::Handle<LLSQLatticeMultiData> > ret; 

      // now loop over all the elements and separate them by lorentz component 
      for(internal_it = internal.begin(); internal_it != internal.end(); ++internal_it)
      {

        // generate the subduced guy with a cartesian index
        generateCartesianRedstarXML tmp(*internal_it); 

        // run any symmetry operations that we might want
        if(m_ini.threePointCorrXMLIni.gParitySymmetry)
          tmp.run_g_parity_symmetry();
        if(m_ini.threePointCorrXMLIni.cubicSymmetry)
          tmp.run_cubic_symmetry();

        // only keep active data
        if(tmp.t.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(0,tmp.t,*internal_it,m_ini));
        if(tmp.x.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(1,tmp.x,*internal_it,m_ini));
        if(tmp.y.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(2,tmp.y,*internal_it,m_ini));
        if(tmp.z.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(3,tmp.z,*internal_it,m_ini));
      }

      // sort them by qsq
      sorted_elements = sortByQ2(unsorted_elements);


#ifdef BUILD_CORRELATORS_PARALLEL
      ret = build_llsq_data_parallel(sorted_elements,m_ini);
#else  

      // loop and get the actual numbers to load into the llsq
      for(sorted_it = sorted_elements.begin(); sorted_it != sorted_elements.end(); ++sorted_it)
      {
        std::cout << "Working on Q2 = " << sorted_it->first << std::endl; 
        std::pair<bool,ADAT::Handle<LLSQLatticeMultiData> > val; 
        val = build_llsq_data(sorted_it->second,m_ini);

        if(val.first)
          ret.push_back(val.second);
        else
          std::cout << "removing Q2 = " << sorted_it->first  << " from llsq (no data) " << std::endl;
      }

#endif

      // dump the bad/missing xml lists
      bad_data_repository.dump_baddies(); 
      //    std::cout << __func__ << ": have " << ret.size() << " different values of Q^2" << std::endl;


      return ret; 
    }



  std::vector<Hadron::KeyHadronNPartNPtCorr_t>
    BuildCorrelators::build_correlator_xml(void)
    {
      POW2_ASSERT(have_ini); 

      // need thes
      registerSubductionTables(); 

      // transform the xml into the internal helicity representation
      std::vector<gParityWorld::GParityHelicityMatrixElement> internal
        = gParityWorld::getGParityHelicityMatrixElementFromXML(
            m_ini.threePointCorrXMLIni.continuumMatElemXML); 

      std::vector<gParityWorld::GParityHelicityMatrixElement>::const_iterator internal_it;
      std::vector<cartesianMatrixElementXML> unsorted_elements; 
      std::vector<cartesianMatrixElementXML>::const_iterator unsorted_iterator; 

      // now loop over all the elements and separate them by lorentz component 
      for(internal_it = internal.begin(); internal_it != internal.end(); ++internal_it)
      {

        // generate the subduced guy with a cartesian index
        generateCartesianRedstarXML tmp(*internal_it); 

        // run any symmetry operations that we might want
        if(m_ini.threePointCorrXMLIni.gParitySymmetry)
          tmp.run_g_parity_symmetry();
        if(m_ini.threePointCorrXMLIni.cubicSymmetry)
          tmp.run_cubic_symmetry();

        // only keep active data
        if(tmp.t.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(0,tmp.t,*internal_it,m_ini));
        if(tmp.x.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(1,tmp.x,*internal_it,m_ini));
        if(tmp.y.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(2,tmp.y,*internal_it,m_ini));
        if(tmp.z.first)
          unsorted_elements.push_back(cartesianMatrixElementXML(3,tmp.z,*internal_it,m_ini));
      }


      // build a map of unique correlators

      std::map<std::string , Hadron::KeyHadronNPartNPtCorr_t > npoint_map; 
      std::map<std::string , Hadron::KeyHadronNPartNPtCorr_t >::const_iterator npoint_it; 
      cartesianMatrixElementXML::listNPointKey::const_iterator cart_it; 

      for(unsorted_iterator = unsorted_elements.begin(); 
          unsorted_iterator != unsorted_elements.end(); 
          ++unsorted_iterator)
        for(cart_it = unsorted_iterator->my_npoint.begin(); 
            cart_it != unsorted_iterator->my_npoint.end(); 
            ++cart_it)
          if (npoint_map.find(Hadron::ensemFileName(cart_it->m_obj)) == npoint_map.end())
            npoint_map[Hadron::ensemFileName(cart_it->m_obj)] = cart_it->m_obj; 

      std::vector<Hadron::KeyHadronNPartNPtCorr_t> ret; 

      for(npoint_it = npoint_map.begin(); npoint_it != npoint_map.end(); ++npoint_it)
        ret.push_back(npoint_it->second); 

      return ret;
    }






} // radmat



#undef BUILD_CORRELATORS_PARALLEL




