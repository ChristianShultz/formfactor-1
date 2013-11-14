/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : construct_correlators.cc

 * Purpose :

 * Creation Date : 04-12-2012

 * Last Modified : Thu 14 Nov 2013 12:19:05 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "construct_correlators.h"
#include "construct_correlators_utils.h"
#include "construct_correlators_bad_data_repository.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "adat/adat_stopwatch.h"
#include "adat/map_obj.h"
#include <string>
#include <vector>
#include <utility>
#include <omp.h>

#define CONSTRUCT_CORRELATORS_PARALLEL
#define TIME_CONSTRUCT_SINGLE_CORRS
#define TIME_CONSTRUCT_ALL_CORRS

namespace radmat
{

  namespace 
  {
    ///////////////////////////////////////////////////////
    BuildCorrsLocalBadDataRepo_t local_bad_data_repo; 


    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    std::pair<bool,rHandle<LLSQLatticeMultiData> >
      build_a_correlator(
          const double qsq,
          const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
          const std::string &sink_id, 
          const std::string &source_id, 
          const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
          const DatabaseInterface_t &db )
      {

#ifdef TIME_CONSTRUCT_SINGLE_CORRS
        Util::StopWatch snoop;
        snoop.start();
#endif

        std::cout << __func__ << ": working on Q2 =" << qsq << std::endl;

        std::pair<bool,std::vector<ConstructCorrsMatrixElement> > tmp; 
        std::vector<Hadron::KeyHadronNPartNPtCorr_t> missed_xml; 
        std::vector<RadmatExtendedKeyHadronNPartIrrep_t> missed_norm;

        tmp = build_correlators(corrs,missed_xml,missed_norm,sink_id,source_id,Z_V,db); 

        if ( !!! missed_xml.empty() ) 
          local_bad_data_repo.insert(omp_get_thread_num(),missed_xml);
        if ( !!! missed_norm.empty() ) 
          local_bad_data_repo.insert(omp_get_thread_num(),missed_norm);

        std::vector<ConstructCorrsMatrixElement>::const_iterator it; 
        bool any_data = false; 
        rHandle<LLSQLatticeMultiData> data(new LLSQLatticeMultiData()); 

        if ( !!! tmp.first ) 
          return std::pair<bool,rHandle<LLSQLatticeMultiData> >(false,data); 

        // another data sanity check.. 
        for(it = tmp.second.begin(); it != tmp.second.end(); ++it)
          if ( it->success ) 
          {
            data->append_row_ensem(it->data,it->tag);  
            any_data = true;
          }

#ifdef TIME_CONSTRUCT_SINGLE_CORRS
        snoop.stop();
        if(any_data)
          std::cout << __func__ << ": q2 = " << qsq 
            << " NxM -> " << data->nrows() << "x" << data->ncols()
            << " took " << snoop.getTimeInSeconds() << " seconds" << std::endl;
#endif

        return std::pair<bool,rHandle<LLSQLatticeMultiData> >(any_data,data); 
      }



    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    std::vector< 
      std::pair<bool,rHandle<LLSQLatticeMultiData> > 
      >
      build_llsq_corrs(
          const std::vector<std::pair<double,std::vector<TaggedEnsemRedstarNPtBlock> > > &loop_data,
          const std::string &snk_id, 
          const std::string &src_id,
          const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
          const DatabaseInterface_t &db)
      {

        std::vector<std::pair<bool,rHandle<LLSQLatticeMultiData> > > ret(loop_data.size()); 
        int idx;
        int sz ( loop_data.size() ) ; 

#ifdef CONSTRUCT_CORRELATORS_PARALLEL
#pragma omp parallel for shared(idx) schedule(dynamic,1)
#endif

        for(idx = 0; idx < sz; ++idx)
          ret[idx] = build_a_correlator( loop_data[idx].first, loop_data[idx].second, 
              snk_id, src_id, Z_V, db);

#ifdef CONSTRUCT_CORRELATORS_PARALLEL
#pragma omp barrier
#endif 

      }



    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // the brain 
    std::vector<rHandle<LLSQLatticeMultiData> >
      do_work(const ThreePointCorrIni_t &ini)
      {

#ifdef TIME_CONSTRUCT_ALL_CORRS
        Util::StopWatch snoop;
        snoop.start(); 
#endif

        double p_factor = mom_factor(ini.xi,ini.L_s); 
        std::string elem_id = ini.matElemID; 
        const radmatDBProp_t *db_prop = &ini.radmatDBProp; 
        const ThreePointCorrXMLIni_t *three_pt = &ini.threePointCorrXMLIni; 

        std::vector<TaggedEnsemRedstarNPtBlock> unsorted_elems; 
        unsorted_elems = tag_lattice_xml(
            &(three_pt->redstar), 
            p_factor, 
            three_pt->maSink, 
            three_pt->maSource, 
            elem_id); 

        std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > sorted_elems;
        std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> >::const_iterator it;
        sorted_elems = sort_tagged_corrs_by_Q2(unsorted_elems); 

        std::vector<std::pair<double,std::vector<TaggedEnsemRedstarNPtBlock> > > loop_data; 

        for(it = sorted_elems.begin(); it != sorted_elems.end(); ++it)
          loop_data.push_back(std::pair<double,std::vector<TaggedEnsemRedstarNPtBlock> >(it->first,it->second)); 


        std::vector<std::pair<bool,rHandle<LLSQLatticeMultiData> > > lattice_data; 
        DatabaseInterface_t db(*db_prop) ; 

        lattice_data = build_llsq_corrs(loop_data,three_pt->sink_id, three_pt->source_id, 
            three_pt->renormalization, db); 


        std::vector<rHandle<LLSQLatticeMultiData> > ret; 
        std::vector<std::pair<bool,rHandle<LLSQLatticeMultiData> > >::const_iterator dcheck; 

        for(dcheck = lattice_data.begin(); dcheck != lattice_data.end(); ++dcheck) 
        {
          if( dcheck->first )
          {
            continue;    
          }
          ret.push_back(dcheck->second); 
        }

        local_bad_data_repo.dump_baddies();



#ifdef TIME_CONSTRUCT_ALL_CORRS
        snoop.stop(); 
        std::cout << " ** time to build all correlators " << snoop.getTimeInSeconds() << " seconds" << std::endl;
#endif
        return ret;  
      }

  } // anonomyous 


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  
  namespace
  {
    void print_timeslice_info(const rHandle<AbsRedstarMergeNPt> &r)
    {
      ADATXML::Array<int>  tslice = r->timeslice_info(); 
      for(int i = 0; i < tslice.size(); ++i)
        std::cout << "npt[" << i << "] lives on timeslice " 
          << tslice[i] << std::endl;  
    }

    void print_npoint_xml_info(const NPointXML & npt)
    {
      std::cout << "NPointXML" << std::endl; 
      std::cout << "version  -> " << npt.version << std::endl; 
      std::cout << "N        -> " << npt.N << std::endl; 
      std::cout << "ensemble -> " << npt.ensemble << std::endl; 
      for(int i =0; i < npt.npoint.size(); ++i)
        std::cout << npt.npoint[i].object_name << ":\n"
          << npt.npoint[i].param->write() << std::endl; 
    }

    void print_merge_npt_data_info(const AbsRedstarMergeNPtData_t & data)
    {
      std::cout << "AbsRedstarMergeNPtData_t " << std::endl; 
      std::cout << "N_EnsemRedstarNPtBlocks = " << data.npoint.size() << std::endl; 
    }
  }

  ///////////////////////////////////////////////////////
  // print some of the xml info for debuggin 
  void
    ConstructCorrelators::print_redstar(const AbstractMergeNamedObject &r)
    {
      rHandle<AbsRedstarMergeNPt> rr = r.param; 
      std::cout << __func__ << ": type -> " << rr->type() << std::endl; 
      print_npoint_xml_info(rr->nptXML());
      print_timeslice_info(rr); 
      print_merge_npt_data_info(rr->data()); 
    }

  ///////////////////////////////////////////////////////
  // construct lots of correlators
  std::vector<rHandle<LLSQLatticeMultiData> >
    ConstructCorrelators::construct_multi_correlators(void) const
    {
      return do_work(m_ini); 
    }


  ///////////////////////////////////////////////////////
  // just make some xml 
  std::vector<Hadron::KeyHadronNPartNPtCorr_t>
    ConstructCorrelators::construct_correlator_xml(void) const
    {

#ifdef TIME_CONSTRUCT_ALL_CORRS
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      double p_factor = mom_factor(m_ini.xi,m_ini.L_s); 
      std::string elem_id = m_ini.matElemID; 
      const radmatDBProp_t *db_prop = &m_ini.radmatDBProp; 
      const ThreePointCorrXMLIni_t *three_pt = &m_ini.threePointCorrXMLIni; 

      std::vector<TaggedEnsemRedstarNPtBlock> unsorted_elems; 
      unsorted_elems = tag_lattice_xml(
          &(three_pt->redstar),
          p_factor, 
          three_pt->maSink, 
          three_pt->maSource, 
          elem_id); 

      std::cout << __PRETTY_FUNCTION__ << ": $#unsorted = " << unsorted_elems.size() << std::endl; 

      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block;
      EnsemRedstarNPtBlock::const_iterator npt; 
      ADAT::MapObject<EnsemRedstarNPtBlock::Obj_t,int> hash; 

      for(block = unsorted_elems.begin(); block != unsorted_elems.end(); ++block)
        for(npt = block->coeff_lattice_xml.begin(); npt != block->coeff_lattice_xml.end(); ++npt)
          if(!!! hash.exist(npt->m_obj) )
            hash.insert(npt->m_obj,1); 

#ifdef TIME_CONSTRUCT_ALL_CORRS
      snoop.stop(); 
      std::cout << " ** time to build xml for " <<  hash.size() 
        << " correlators " << snoop.getTimeInSeconds() << " seconds" << std::endl;
#endif

      return hash.keys(); 
    }

} // radmat


