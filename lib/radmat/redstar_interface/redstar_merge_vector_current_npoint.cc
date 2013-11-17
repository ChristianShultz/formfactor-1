/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_merge_vector_current_npoint.cc

 * Purpose :

 * Creation Date : 12-11-2013

 * Last Modified : Fri 15 Nov 2013 06:34:44 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_merge_vector_current_npoint.h"
#include "redstar_abstract_merge_npoint.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/pow2assert.h"
#include "redstar_single_particle_meson_block.h"
#include "redstar_unimproved_vector_current.h"
#include "adat/adat_stopwatch.h"
#include "hadron/irrep_util.h"
#include "hadron/ensem_filenames.h"


#define  DEBUG_MSG_OFF
#define  DEBUG_HANDLE_ON
#include "debug_props.h"

namespace
{
  int SINK_INDEX(0);
  int INSERTION_INDEX(1);
  int SOURCE_INDEX(2);
  int SINK_INDEX_D1(1);
  int INSERTION_INDEX_D1(2);
  int SOURCE_INDEX_D1(3);
}

namespace radmat
{

  namespace 
  {

    ////////////////////////////////////////////////////
    // some sanity checks  
    void check_n(const int N)
    {
      POW2_ASSERT(N == 3); 
    }

    ////////////////////////////////////////////////////
    void check_version(const int v)
    {
      POW2_ASSERT(v == 0); 
    }

    ////////////////////////////////////////////////////
    void check_meson(const rHandle<AbsRedstarBlock_t> p)
    {
      bool meson = false; 

      meson |= (p->type() == Stringify<RedstarSingleParticleMesonBlock>()); 

      POW2_ASSERT( meson ); 
    }

    ////////////////////////////////////////////////////
    void check_photon(const rHandle<AbsRedstarBlock_t> p)
    {
      bool photon = false; 

      photon |= (p->type() == Stringify<RedstarUnimprovedVectorCurrentBlock>());

      POW2_ASSERT( photon ); 
    }

    ////////////////////////////////////////////////////
    // do the cut
    bool cutMomentum(const ADATXML::Array<int> &mom,
        const int minmom ,
        const int maxmom)
    {
      int sq(0);
      POW2_ASSERT(mom.size() == 3); 
      for(int i = 0; i < 3; ++i)
        sq += mom[i]*mom[i];


      //      std::cout << __func__ << ": " << sq 
      //        << " " << minmom << " " << maxmom 
      //        << " " << ( (sq >= minmom) && (sq <= maxmom) ) << std::endl;

      return ( (sq >= minmom) && (sq <= maxmom) ) ;  // false if its too big or too small 
    }

    ////////////////////////////////////////////////////
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

    ////////////////////////////////////////////////////
    // determine the momentum transfer -- this is our phase convention, 
    //    mess with it at your own risk
    ADATXML::Array<int> momentumTransfer(const ADATXML::Array<int> &sink, 
        const ADATXML::Array<int> &source, 
        const bool create)
    { 
      POW2_ASSERT(sink.size() == source.size());
      POW2_ASSERT(sink.size() == 3); 
      ADATXML::Array<int> ret;
      ret.resize(3);

      ret[0] = source[0] - sink[0];
      ret[1] = source[1] - sink[1];
      ret[2] = source[2] - sink[2];

      if (create) 
        return -ret;

      std::cout << "pf" << sink[0] << sink[1] << sink[2] 
        << "  pi" << source[0] << source[1] << source[2] 
        << "   q" << ret[0] << ret[1] << ret[2]; 

      return ret; 
    }

    ////////////////////////////////////////////////////
    // upcast to find data
    ADATXML::Array<int> 
      find_meson_momentum(const rHandle<AbsRedstarInput_t> foo)
      {
        DEBUG_MSG(entering); 
        DEBUG_HANDLE(foo); 

        ADATXML::Array<int> ret; 
        if(foo->type() == Stringify<RedstarSingleParticleMesonInput>())
        {
          DEBUG_MSG(casting);
          rHandle<RedstarSingleParticleMesonInput> p(foo); 
          DEBUG_HANDLE(p);
          ret = p->mom; 
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
            << ": Error, unrecognized meson " << foo->type() << std::endl; 
          exit(1); 
        }

        DEBUG_MSG(exiting);

        return ret; 
      }

    ////////////////////////////////////////////////////
    // upcast to find data
    ADATXML::Array<int> 
      find_photon_momentum(const rHandle<AbsRedstarInput_t> f)
      {
        DEBUG_MSG(entering); 
        DEBUG_HANDLE(f); 

        ADATXML::Array<int> ret; 
        if ( f->type() == Stringify<RedstarUnimprovedVectorCurrentInput>())
        {
          DEBUG_MSG(casting); 
          rHandle<RedstarUnimprovedVectorCurrentInput> p(f); 
          DEBUG_HANDLE(p);
          ret = p->mom;     
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
            << ": Error, unrecognized photon " << f->type() << std::endl; 
          exit(1); 
        }

        DEBUG_MSG(exiting);

        return ret; 
      }


    ////////////////////////////////////////////////////
    // upcast to set data
    void set_q(rHandle<AbsRedstarInput_t> f, 
        const ADATXML::Array<int> psnk,
        const ADATXML::Array<int> psrc)
    {
      DEBUG_MSG(entering);
      if ( f->type() == Stringify<RedstarUnimprovedVectorCurrentInput>())
      {
        DEBUG_MSG(casting); 
        rHandle<RedstarUnimprovedVectorCurrentInput> p(f); 
        DEBUG_HANDLE(p);
        p->mom = momentumTransfer(psnk,psrc,p->creation_op);     
      }
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
          << ": Error, unrecognized insertion " << f->type() << std::endl; 
        exit(1); 
      }
      DEBUG_MSG(exiting);
    }


    ////////////////////////////////////////////////////
    //  conserve momentum  
    void fill_insertion_momentum(std::vector<rHandle<AbsRedstarInput_t> > &inp)
    {
      DEBUG_MSG(entering);

      ADATXML::Array<int> psnk, psrc,q;  

      psnk = find_meson_momentum(inp[SINK_INDEX]); 
      psrc = find_meson_momentum(inp[SOURCE_INDEX]); 

      set_q(inp[INSERTION_INDEX],psnk,psrc);  

      DEBUG_MSG(exiting);
    } 



    ////////////////////////////////////////////////////
    // hardwire for three point
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
        npt.npoint[SINK_INDEX_D1] = snk.m_obj;
        npt.npoint[INSERTION_INDEX_D1] = ins.m_obj; 
        npt.npoint[SOURCE_INDEX_D1] = src.m_obj;
        npt.ensemble = ensemble;

        ret.m_coeff = coeff; 
        ret.m_obj = npt; 
        return ret; 
      }



    ////////////////////////////////////////////////////
    // sum this thing out 
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


    ////////////////////////////////////////////////////
    // 0.7 * A + 0.5 *B + 0.2*A -> 0.9*A + 0.5*B
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

    ////////////////////////////////////////////////////
    // an intermediary 
    struct merge
    {
      merge(void) : success(false) {}

      bool success; 
      EnsemRedstarNPtBlock npt; 
      std::vector<rHandle<AbsRedstarInput_t> > input; 
    };

    ////////////////////////////////////////////////////
    // a true means continue working with this combo
    bool allowed_photon_momentum(const rHandle<AbsRedstarInput_t> photon , 
        const AbstractBlockNamedObject &ins)
    {

      bool success = true; 
      ADATXML::Array<int> q; 
      int pmin,pmax;  

      if(ins.param->type() == Stringify<RedstarUnimprovedVectorCurrentXML>())
      {
        rHandle<RedstarUnimprovedVectorCurrentXML> p(ins.param); 
        DEBUG_HANDLE(p);
        pmin = p->pmin; 
        pmax = p->pmax; 

        q = find_photon_momentum(photon) ;  
      }
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << __FILE__ 
          << ": Error, unrecognized insertionXML " << ins.param->type() << std::endl; 
        exit(1); 
      }

      success &= cutMomentum(q,pmin,pmax);
      // need to put a check here or else the 
      // irrep utils will bug out on weird things
      // like 2 -2 1 
      if(success)  
        success &= acceptable_mom(q); 

      return success; 
    }

    ////////////////////////////////////////////////////
    // the work horse -- handle subducing  
    merge do_merge(
        const rHandle<AbsRedstarBlock_t> snk, 
        const rHandle<AbsRedstarInput_t> snki,
        const rHandle<AbsRedstarBlock_t> ins, 
        const rHandle<AbsRedstarInput_t> insi,
        const rHandle<AbsRedstarBlock_t> src, 
        const rHandle<AbsRedstarInput_t> srci,
        const AbstractBlockNamedObject &insXML,
        const std::string ensemble) 
    {
      DEBUG_MSG(entering);

      merge ret; 
      EnsemRedstarNPtBlock npt; 
      std::vector< rHandle<AbsRedstarInput_t > > input(3); 
      bool success(true); 

      DEBUG_HANDLE(snk); 
      DEBUG_HANDLE(snki); 
      DEBUG_HANDLE(ins); 
      DEBUG_HANDLE(insi); 
      DEBUG_HANDLE(src); 
      DEBUG_HANDLE(srci); 

      // make unique copies of the input data and shove it into 
      // THIS vector which corresponds to the xml for THIS mat
      // elem, then conserve momentum, all further operations 
      // corresponding to this mat elem must use the input vector
      // else we end up using uninitialized variables which sucks
      input[SINK_INDEX] = rHandle<AbsRedstarInput_t>(snki->clone()); 
      input[INSERTION_INDEX] = rHandle<AbsRedstarInput_t>(insi->clone()); 
      input[SOURCE_INDEX] = rHandle<AbsRedstarInput_t>(srci->clone()); 

      // conserve momentum with the photon
      fill_insertion_momentum(input); 

      // check that it is an allowed photon momentum 
      success &= allowed_photon_momentum( input[INSERTION_INDEX] , insXML ); 

      if ( success ) 
      {
        DEBUG_MSG(subducing);

        // do the subduction -- need to use input vector
        //    since he has had the photon momentum filled 
        EnsemRedstarBlock snkb,insb,srcb; 
        snkb = (*snk)(input[SINK_INDEX]); 
        insb = (*ins)(input[INSERTION_INDEX]); 
        srcb = (*src)(input[SOURCE_INDEX]); 

        npt = sum_duplicates( loop_npt(snkb,insb,srcb,ensemble) ) ; 
        success &= (npt.size() != 0); 
      }

      ret.npt = npt; 
      ret.input = input; 
      ret.success = success ; 

      DEBUG_MSG(exiting);
      return ret; 
    }



    ////////////////////////////////////////////////////
    // try to build ALL xml data
    AbsRedstarMergeNPtData_t handle_work(const NPointXML &npt)
    {
      DEBUG_MSG(entering); 

      // A RGE-ism 
      Util::StopWatch snoop; 
      snoop.start(); 

      std::vector<EnsemRedstarNPtBlock> npt_xml;
      std::vector<std::vector< rHandle< AbsRedstarInput_t > > > npt_inp; 
      AbsRedstarMergeNPtData_t ret; 

      // is it a three point and do we understand what it is
      check_n(npt.N); 
      check_version(npt.version); 

      // iteration variables, pull out some named data
      AbsRedstarXMLInterface_t::const_iterator src,ins,snk;
      AbstractBlockNamedObject source,insertion,sink; 
      sink = npt.npoint[SINK_INDEX];
      insertion = npt.npoint[INSERTION_INDEX]; 
      source = npt.npoint[SOURCE_INDEX];

#ifdef DEBUG_MSG_ON 
      std::cout << __PRETTY_FUNCTION__ << " :: type_info \n" 
        << "sink -> " << sink.param->type() 
        << "\ninsertion -> " << insertion.param->type() 
        << "\nsource -> " << source.param->type() << std::endl;
#endif

      // pointers to the functors we need to use
      rHandle<AbsRedstarBlock_t> snkF = sink.param->objFunctorPtr;
      rHandle<AbsRedstarBlock_t> insF = insertion.param->objFunctorPtr;
      rHandle<AbsRedstarBlock_t> srcF = source.param->objFunctorPtr;

      DEBUG_HANDLE(snkF);
      DEBUG_HANDLE(insF); 
      DEBUG_HANDLE(srcF); 


      // check that they do correspond to mesons/photon 
      check_meson(snkF); 
      check_meson(srcF); 
      check_photon(insF); 

      DEBUG_MSG(bigloop);

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
                insertion,
                npt.ensemble);

            // chuck conditions
            if ( !!! tmp.success ) 
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

      DEBUG_MSG(exiting);

      return ret;
    }

  } // anonomyous 


  ////////////////////////////////////////////////////
  // wrapper
  void RedstarMergeVectorCurrentThreePoint::do_work(void)  
  {
    my_data = handle_work(my_npt); 
  }

} // radmat



