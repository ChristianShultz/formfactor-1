/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : construct_correlators_utils.cc

 * Purpose :

 * Creation Date : 13-11-2013

 * Last Modified : Wed 13 Nov 2013 05:24:55 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "construct_correlators_utils.h"
#include "radmat/fake_data/fake_3pt_function_aux.h"
#include "semble/semble_meta.h"
#include "semble/semble_file_management.h"
#include "adat/adat_stopwatch.h"
#include "hadron/ensem_filenames.h"

#define DO_TIMING_SORT_MAT_ELEMS_BY_Q2
#define DO_TIMING_CACHE_NORM_MAT_ELEMS
#define DO_TIMING_SUM_NORM_MAT_ELEMS

namespace radmat
{


  namespace
  {

    typedef ADAT::MapObject<Hadron::KeyHadronNPartNPtCorr_t,
            ENSEM::EnsemVectorComplex> singleThreadQ2NormalizedCorrCache;

  ////////////////////////////////////////////////////////////////////

    // do the work and return the unique elems from the vector
    singleThreadQ2NormalizedCorrCache
      normalize_correlators( 
          const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
          std::vector<Hadron::KeyHadronNPartNPtCorr_t> &missed_xml, 
          std::vector<RadmatExtendedKeyHadronNPartIrrep_t> &missed_norm, 
          const std::string &sink_id, 
          const std::string &source_id, 
          const ThreePointCorrXMLIni_t::RenormalizationProp & Z_V,
          const DatabaseInterface_t &db
          )
      {
#ifdef DO_TIMING_CACHE_NORM_MAT_ELEMS
        Util::StopWatch snoop; 
        snoop.start();
#endif

        singleThreadQ2NormalizedCorrCache cache; 
        std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block; 
        EnsemRedstarNPtBlock::const_iterator npt; 

        // set up a zero
        ENSEM::EnsemVectorComplex zero; 
        zero.resize(1);
        zero.resizeObs(1);
        zero = SEMBLE::toScalar(0.); 

        for ( block = corrs.begin(); block != corrs.end(); ++block)
          for(npt = block->coeff_lattice_xml.begin(); npt != block->coeff_lattice_xml.end(); ++npt) 
          {
            if(cache.exist(npt->m_obj))
              continue; 

            // can we make it?
            bool found = true; 
            DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
            if( !!!db.exists(npt->m_obj) ) 
            {
              missed_xml.push_back(npt->m_obj); 
              found = false;   
            }
            if( !!!db.exists(k.sink) ) 
            {
              missed_norm.push_back(k.sink); 
              found = false; 
            }
            if( !!!db.exists(k.source) ) 
            {
              missed_norm.push_back(k.source); 
              found = false; 
            }

            // break early if we're missing ingredients 
            if ( !!! found ) 
              continue; 

            // --- DO NORMALIZATION/RENORMALIZATION HERE
            // these are real but the * operator is only overloaded
            // for complex types so shove a zero down its throat
            // to get the desired result
            double Z = Z_V.RGE_prop;  
            ENSEM::Complex Z_c; 

            if ( block->continuum_tag.mu() == 4 ) 
            {
              Z_c = SEMBLE::toScalar(std::complex<double>(Z*Z_V.Z_t)); 
            }
            else if ( (block->continuum_tag.mu() > 0)
                && (block->continuum_tag.mu() < 4) )
            {
              Z_c = SEMBLE::toScalar(std::complex<double>(Z*Z_V.Z_s)); 
            }
            else
            {
              std::cerr << __PRETTY_FUNCTION__ 
                << ": Error, I don't know what lorentz index " 
                << block->continuum_tag.mu() << " means " 
                << std::endl; 
              exit(1);  
            }

            // initialize some variables        
            ThreePtPropagationFactor<double> propagation_factor;

            ENSEM::EnsemVectorComplex corr_tmp = db.fetch(npt->m_obj);
            RadmatMassOverlapData_t source = db.fetch(k.source); 
            RadmatMassOverlapData_t sink = db.fetch(k.sink); 

            std::string outstem = Hadron::ensemFileName(npt->m_obj); 
            std::string pth = SEMBLE::SEMBLEIO::getPath();
            std::stringstream path;

            path << pth << "Q2_" << block->qsq_tag(); 
            SEMBLE::SEMBLEIO::makeDirectoryPath(path.str()); 

            path << "/correlator_normalization";
            SEMBLE::SEMBLEIO::makeDirectoryPath(path.str());

            path << "/" << outstem;
            ENSEM::write(path.str() + std::string("_corr_pre") , corr_tmp); 
            ENSEM::EnsemVectorComplex norm = corr_tmp * SEMBLE::toScalar(0.); 

            // the hadron key uses 1 based arrays
            // NB: assumption that npt is organized like <sink, ins , source>
            const int t_sink(npt->m_obj.npoint[1].t_slice); 
            const int t_source(npt->m_obj.npoint[3].t_slice);

            // sanity
            POW2_ASSERT(t_source < t_sink); 

            // NB: the indexing here assumes [tsource,tsink] ie: inclusive range
            for(int t_ins = t_source; t_ins <= t_sink; ++t_ins)
            {
              ENSEM::EnsemReal prop = propagation_factor(sink.E(),sink.Z(),t_sink,t_ins,
                  source.E(),source.Z(),t_source);

              ENSEM::pokeObs(corr_tmp,ENSEM::peekObs(corr_tmp,t_ins)/prop/Z_c,t_ins);

              ENSEM::pokeObs(norm,prop*Z_c,t_ins); 
            } // end loop over t_ins

            ENSEM::write(path.str() + std::string("_corr_post"), corr_tmp);
            ENSEM::write(path.str() + std::string("_norm") , norm);

            cache.insert(npt->m_obj,corr_tmp); 

          } // end loop over npt  

#ifdef DO_TIMING_CACHE_NORM_MAT_ELEMS
        snoop.stop(); 
        std::cout << __PRETTY_FUNCTION__ << ": normalizing "
          << cache.size() << " unique npts took"
          << snoop.getTimeInSeconds() << " seconds" << std::endl; 
#endif

        return cache; 
      } 

  } // anonomyous 


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  // brute force sort on qsq
  std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > 
    sort_tagged_corrs_by_Q2(const std::vector<TaggedEnsemRedstarNPtBlock> &unsorted)
    {

#ifdef DO_TIMING_SORT_MAT_ELEMS_BY_Q2
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> > ret; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator it; 

      // loop everything and toss it into a vector if it matches else make a new vector
      for(it = unsorted.begin(); it != unsorted.end(); ++it)
      {

        if(ret.find(it->qsq_tag()) != ret.end())
          ret.find(it->qsq_tag())->second.push_back(*it);
        else
          ret.insert(std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> >::value_type(
                it->qsq_tag(),std::vector<TaggedEnsemRedstarNPtBlock>(1,*it))); 
      }


#ifdef DO_TIMING_SORT_MAT_ELEMS_BY_Q2
      snoop.stop(); 

      std::map<double,std::vector<TaggedEnsemRedstarNPtBlock> >::const_iterator it1; 
      int ct(0);
      for(it1 = ret.begin(); it1 != ret.end(); ++it1)
        ct += it1->second.size();

      std::cout << __PRETTY_FUNCTION__ << ": sorting, " << ct << " elems took " 
        << snoop.getTimeInSeconds() << " seconds "  << std::endl; 
#endif 
      return ret;
    }


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  //build a correlator
  std::pair<bool, std::vector<ConstructCorrsMatrixElement> > 
    build_correlators(
        const std::vector<TaggedEnsemRedstarNPtBlock> &corrs,
        std::vector<Hadron::KeyHadronNPartNPtCorr_t> &missed_xml, 
        std::vector<RadmatExtendedKeyHadronNPartIrrep_t> &missed_norm, 
        const std::string &sink_id, 
        const std::string &source_id, 
        const ThreePointCorrXMLIni_t::RenormalizationProp &Z_V,
        const DatabaseInterface_t &db)
    {
#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      Util::StopWatch snoop;
      snoop.start(); 
#endif

      std::vector<ConstructCorrsMatrixElement> ret; 
      double qsq; 

      singleThreadQ2NormalizedCorrCache npoint_cache; 
      npoint_cache = normalize_correlators( corrs, missed_xml, missed_norm, 
          sink_id, source_id, Z_V, db);

      ENSEM::EnsemReal ScalarZeroR;
      ENSEM::EnsemVectorComplex VectorZeroC; 
      bool any_data = false; 
      std::vector<TaggedEnsemRedstarNPtBlock>::const_iterator block; 
      EnsemRedstarNPtBlock::const_iterator npt; 

      bool found = false; 
      for ( block = corrs.begin(); block != corrs.end(); ++block)
      {
        if (found ) 
          break; 

        for(npt = block->coeff_lattice_xml.begin(); npt != block->coeff_lattice_xml.end(); ++npt) 
          if(npoint_cache.exist(npt->m_obj))
          {
            found = true; 
            DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
            RadmatMassOverlapData_t dummy = db.fetch(k.source); 
            ScalarZeroR = SEMBLE::toScalar(0.) * dummy.E();  
            VectorZeroC = SEMBLE::toScalar(std::complex<double>(0.,0.))* db.fetch(npt->m_obj); 
            break; 
          }
      }     

      // can't do anything -- no data 
      if ( !!! found ) 
        return std::pair<bool, std::vector<ConstructCorrsMatrixElement> >(false,ret);  

      for ( block = corrs.begin(); block != corrs.end(); ++block)
      {
        ENSEM::EnsemVectorComplex corr; 
        ENSEM::EnsemReal E_snk, E_src; 

        corr = VectorZeroC; 
        E_snk = ScalarZeroR; 
        E_src = ScalarZeroR;
        int ct(0); 
        bool success = true;  

        for(npt = block->coeff_lattice_xml.begin(); npt != block->coeff_lattice_xml.end(); ++npt) 
        {
          success &= npoint_cache.exist(npt->m_obj);

          // abandon rest of loop if we miss one
          if ( !!! success ) 
            break; 

          DatabaseInterface_k k(npt->m_obj,sink_id,source_id);  
          RadmatMassOverlapData_t source = db.fetch(k.source); 
          RadmatMassOverlapData_t sink = db.fetch(k.sink); 

          E_snk = E_snk + sink.E(); 
          E_src = E_src + source.E(); 

          corr = corr + npt->m_coeff * npoint_cache[npt->m_obj];

          ++ct; 
        }

        if ( ct != 0)
        {
          E_snk = E_snk / SEMBLE::toScalar(double(ct)); 
          E_src = E_src / SEMBLE::toScalar(double(ct)); 
        }

        LatticeMultiDataTag tag = block->continuum_tag;
        tag.E_f = E_snk; 
        tag.E_i = E_src; 

        any_data |= success; 

        ret.push_back( ConstructCorrsMatrixElement(corr,tag,success) );
      }

#ifdef DO_TIMING_SUM_NORM_MAT_ELEMS
      snoop.stop(); 
      std::cout << __PRETTY_FUNCTION__ << " :qsq = " << qsq
        << " :building " << ret.size() << " elems took "
        << snoop.getTimeInSeconds() << " seconds" << std::endl; 
#endif

      return std::pair<bool, std::vector<ConstructCorrsMatrixElement> >(any_data,ret);  
    }


} // radmat 

