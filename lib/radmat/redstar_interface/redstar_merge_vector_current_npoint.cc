/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_merge_vector_current_npoint.cc

 * Purpose :

 * Creation Date : 12-11-2013

 * Last Modified : Thu 14 Nov 2013 10:34:14 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_merge_vector_current_npoint.h"
#include "redstar_abstract_merge_npoint.h"
#include "handle_interface.h"
#include "radmat/utils/pow2assert.h"
#include "redstar_single_particle_meson_block.h"
#include "redstar_unimproved_vector_current.h"
#include "adat/adat_stopwatch.h"
#include "hadron/irrep_util.h"
#include "hadron/ensem_filenames.h"

namespace radmat
{

  namespace 
  {

    // some sanity checks  
    void check_n(const int N)
    {
      POW2_ASSERT(N == 3); 
    }

    void check_version(const int v)
    {
      POW2_ASSERT(v == 0); 
    }

    void check_meson(const ADAT::Handle<AbsRedstarBlock_t> p)
    {
      bool meson = false; 

      meson |= (p->type() == Stringify<RedstarSingleParticleMesonBlock>()); 

      POW2_ASSERT( meson ); 
    }

    void check_photon(const ADAT::Handle<AbsRedstarBlock_t> p)
    {
      bool photon = false; 

      photon |= (p->type() == Stringify<RedstarUnimprovedVectorCurrentBlock>());

      POW2_ASSERT( photon ); 
    }


    // do the cut
    bool cutMomentum(const ADATXML::Array<int> &mom, const int minmom , const int maxmom)
    {
      int sq(0);
      for(int i = 0; i < 3; ++i)
        sq += mom[i]*mom[i];
      return   !!! ( (sq >= minmom) && (sq <= maxmom) ) ;  // false if its too big or too small 
    }

    // LG that we know about
    bool acceptable_mom(const ADATXML::Array<int> &mom)
    {
      std::string lg = Hadron::generateLittleGroup(mom); 
      if(lg == "Oh")
        return true;
      if(lg == "D4")
        return true; 
      if(lg == "D2")
        return true; 
      if(lg == "D3")
        return true; 

      std::cout << __func__ << ": not supporting LG " << lg << std::endl;
      return false;
    }

    // determine the momentum transfer -- this is our phase convention, 
    //    mess with it at your own risk
    ADATXML::Array<int> momentumTransfer(const ADATXML::Array<int> &sink, 
        const ADATXML::Array<int> &source, 
        const bool create)
    { 
      ADATXML::Array<int> ret;
      ret.resize(3);

      ret[0] = source[0] - sink[0];
      ret[1] = source[1] - sink[1];
      ret[2] = source[2] - sink[2];

      if (create) 
        return -ret;

      return ret; 
    }

    // upcast to find data
    ADATXML::Array<int> 
      find_meson_momentum(const ADAT::Handle<AbsRedstarInput_t> foo)
      {
        ADATXML::Array<int> ret; 
        if(foo->type() == Stringify<RedstarSingleParticleMesonInput>())
        {
          const RedstarSingleParticleMesonInput *p; 
          p = dynamic_cast_handle<RedstarSingleParticleMesonInput,AbsRedstarInput_t>(foo); 
          ret = p->mom; 
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
            << ": Error, unrecognized meson " << foo->type() << std::endl; 
          exit(1); 
        }
        return ret; 
      }


    // upcast to set data
    void set_q(ADAT::Handle<AbsRedstarInput_t> f, 
        const ADATXML::Array<int> psnk,
        const ADATXML::Array<int> psrc)
    {
      if ( f->type() == Stringify<RedstarUnimprovedVectorCurrentInput>())
      {
        RedstarUnimprovedVectorCurrentInput *p; 
        p = dynamic_cast<RedstarUnimprovedVectorCurrentInput*>(&*f); 
        p->mom = momentumTransfer(psnk,psrc,p->creation_op);     
      }
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
          << ": Error, unrecognized insertion " << f->type() << std::endl; 
        exit(1); 
      }
    }


    //  conserve momentum  
    void fill_insertion_momentum(std::vector<ADAT::Handle<AbsRedstarInput_t> > &inp)
    {
      bool success = true;  
      ADATXML::Array<int> psnk, psrc,q;  

      psnk = find_meson_momentum(inp[0]); 
      psrc = find_meson_momentum(inp[2]); 

      set_q(inp[0],psnk,psrc);  
    } 



    EnsemRedstarNPtBlock::ListObj_t 
      construct_npoint( const EnsemRedstarBlock::ListObj_t snk, 
          const EnsemRedstarBlock::ListObj_t ins, 
          const EnsemRedstarBlock::ListObj_t src,
          const std::string &ensemble) 
      {
        EnsemRedstarNPtBlock::ListObj_t ret; 
        Hadron::KeyHadronNPartNPtCorr_t npt; 
        ENSEM::Complex coeff;        

        coeff = snk.m_coeff * ins.m_coeff * src.m_coeff; 
        npt.npoint.resize(3); 
        npt.npoint[1] = snk.m_obj;
        npt.npoint[2] = ins.m_obj; 
        npt.npoint[3] = src.m_obj;

        ret.m_coeff = coeff; 
        ret.m_obj = npt; 
        return ret; 
      }



    EnsemRedstarNPtBlock loop_npt(const EnsemRedstarBlock &snkl, 
        const EnsemRedstarBlock &insl, 
        const EnsemRedstarBlock &srcl,
        const std::string &ensemble)
    {
      EnsemRedstarNPtBlock ret; 
      EnsemRedstarBlock::const_iterator snk,ins,src; 

      for(snk = snkl.begin(); snk != snkl.end(); ++snk)
        for(ins = insl.begin(); ins != insl.end(); ++ins)
          for(src = srcl.begin(); src != srcl.end(); ++src)
          {
            ret = ret + construct_npoint(*snk,*ins,*src,ensemble); 
          }

      return ret; 
    }


    EnsemRedstarNPtBlock sum_duplicates( const EnsemRedstarNPtBlock &in)
    {
      EnsemRedstarNPtBlock ret; 

      // loop and combine coefficients for duplicate npoints
      std::map<std::string,
        std::pair<EnsemRedstarNPtBlock::Coeff_t,EnsemRedstarNPtBlock::Obj_t> > coeff_map;

      EnsemRedstarNPtBlock::const_iterator it; 
      for(it = in.begin(); it != in.end(); ++it)
      {
        ENSEM::Complex coeff = it->m_coeff;
        std::string key = ensemFileName(it->m_obj); 

        if(coeff_map.find(key) != coeff_map.end())
          coeff += coeff_map[key].first;

        // reinsert
        coeff_map[key] = 
          std::pair<EnsemRedstarNPtBlock::Coeff_t,EnsemRedstarNPtBlock::Obj_t>(coeff,it->m_obj); 
      }

      // write the output list
      std::map<std::string,
        std::pair<EnsemRedstarNPtBlock::Coeff_t,EnsemRedstarNPtBlock::Obj_t>
          >::const_iterator map_it;

      // scan it for zeros in the process
      for(map_it = coeff_map.begin(); map_it != coeff_map.end(); ++map_it)
        if( ENSEM::toDouble(ENSEM::localNorm2(map_it->second.first) ) > 0.0001 )
          ret = ret + EnsemRedstarNPtBlock::ListObj_t(map_it->second.first,map_it->second.second);  

      return ret; 
    }

    // an intermediary 
    struct merge
    {
      EnsemRedstarNPtBlock npt; 
      std::vector<ADAT::Handle<AbsRedstarInput_t> > input; 
    };

    // the work horse 
    merge do_merge(
        const ADAT::Handle<AbsRedstarBlock_t> snk, 
        const ADAT::Handle<AbsRedstarInput_t> snki,
        const ADAT::Handle<AbsRedstarBlock_t> ins, 
        const ADAT::Handle<AbsRedstarInput_t> insi,
        const ADAT::Handle<AbsRedstarBlock_t> src, 
        const ADAT::Handle<AbsRedstarInput_t> srci,
        const std::string ensemble) 
    {
      merge ret; 
      EnsemRedstarNPtBlock npt; 
      std::vector< ADAT::Handle<AbsRedstarInput_t > > input(3); 

      input[0] = snki; 
      input[0] = insi; 
      input[0] = srci; 

      // conserve momentum with the photon
      fill_insertion_momentum(input); 

      // do the subduction 
      EnsemRedstarBlock snkb,insb,srcb; 
      snkb = (*snk)(snki); 
      insb = (*ins)(insi); 
      srcb = (*src)(srci); 

      npt = sum_duplicates( loop_npt(snkb,insb,srcb,ensemble) ) ; 

      ret.npt = npt; 
      ret.input = input; 
      return ret; 
    }


    bool check_merge(const merge &tmp, const AbstractBlockNamedObject &ins)
    {
      if( tmp.npt.size() == 0) 
        return false; 

      bool success = true; 
      ADATXML::Array<int> q; 
      int pmin,pmax;  

      if(ins.param->type() == Stringify<RedstarUnimprovedVectorCurrentXML>())
      {
        const RedstarUnimprovedVectorCurrentXML *p; 
        p = dynamic_cast_handle<RedstarUnimprovedVectorCurrentXML,AbsRedstarXMLInterface_t>(ins.param); 
        pmin = p->pmin; 
        pmax = p->pmax; 

        EnsemRedstarNPtBlock::ListObj_t dum =  *(tmp.npt.begin()); 
        q = dum.m_obj.npoint[2].irrep.mom;  
      }
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
          << ": Error, unrecognized insertionXML " << ins.param->type() << std::endl; 
        exit(1); 
      }

      success &= cutMomentum(q,pmin,pmax);
      success &= acceptable_mom(q); 

      return success; 
    }


    AbsRedstarMergeNPtData_t handle_work(const NPointXML &npt)
    {
      // An RGE-ism 
      Util::StopWatch snoop; 
      snoop.start(); 

      std::vector<EnsemRedstarNPtBlock> npt_xml;
      std::vector<std::vector< ADAT::Handle< AbsRedstarInput_t > > > npt_inp; 
      AbsRedstarMergeNPtData_t ret; 

      // is it a three point and do we understand what it is
      check_n(npt.N); 
      check_version(npt.version); 

      // iteration variables, pull out some named data
      AbsRedstarXMLInterface_t::const_iterator src,ins,snk;
      AbstractBlockNamedObject source,insertion,sink; 
      sink = npt.npoint[0];
      insertion = npt.npoint[1]; 
      source = npt.npoint[2];

      // pointers to the functors we need to use
      ADAT::Handle<AbsRedstarBlock_t> snkF = sink.param->objFunctorPtr;
      ADAT::Handle<AbsRedstarBlock_t> insF = insertion.param->objFunctorPtr;
      ADAT::Handle<AbsRedstarBlock_t> srcF = source.param->objFunctorPtr;

      // check that they do correspond to mesons/photon 
      check_meson(snkF); 
      check_meson(srcF); 
      check_photon(insF); 

      // big loop, each continuum operator is a sum of lattice operators, 
      //    this sum needs to be multiplied out, do it here
      for(snk = sink.param->begin(); snk != sink.param->end(); ++snk)
        for(ins = insertion.param->begin(); ins != insertion.param->end(); ++ins)
          for(src = source.param->begin(); src != source.param->end(); ++src)
          {
            // this continuum combination -- this is the workhorse  
            merge tmp = do_merge(snkF,*snk,
                insF,*ins,
                srcF,*src,
                npt.ensemble);

            // more chuck conditions
            if ( ! check_merge(tmp,insertion) ) 
              continue; 

            // push it to the big lists
            npt_xml.push_back(tmp.npt); 
            npt_inp.push_back(tmp.input); 
          }

      // combine to make indexed data
      ret.npoint = npt_xml; 
      ret.input = npt_inp; 

      // finished 
      snoop.stop(); 
      std::cout << __PRETTY_FUNCTION__ << ":" << __FILE__ 
        << " took " << snoop.getTimeInSeconds() << " seconds" << std::endl; 

      return ret;
    }

  } // anonomyous 


  void RedstarMergeVectorCurrentThreePoint::do_work(void)  
  {
    my_data = handle_work(my_npt); 
  }

} // radmat


