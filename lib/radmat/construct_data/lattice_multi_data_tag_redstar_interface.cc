/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lattice_multi_data_tag_redstar_interface.cc

 * Purpose :

 * Creation Date : 12-11-2013

 * Last Modified : Thu 14 Nov 2013 05:06:35 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "lattice_multi_data_tag_redstar_interface.h"
#include "radmat/utils/handle.h"
#include "radmat/utils/mink_qsq.h"
#include "radmat/utils/pow2assert.h"
#include <sstream>
#include <omp.h>

#define PARALLEL_TAG_REDSTAR_DATA

namespace radmat
{

  namespace 
  {

    // intermediary data storage
    struct 
      data_store
      {
        std::string file_id;
        int jmu; 
        std::string mat_elem_id; 
        ADATXML::Array<int> p_f;
        ADATXML::Array<int> p_i; 
      };

    // actually makes the tag
    LatticeMultiDataTag
      generate_tag(const data_store &d,
          const double mom_factor, 
          const double msnk,
          const double msrc)
      {
        LatticeMultiDataTag ret; 

        ret.file_id = d.file_id; 
        ret.jmu = d.jmu; 
        ret.mat_elem_id = d.mat_elem_id; 
        ret.p_f = d.p_f;
        ret.p_i = d.p_i; 
        ret.mom_fac = mom_factor; 
        ret.set_qsq_label(
            Mink_qsq(ret.p_f,msnk,ret.p_i,msrc,ret.mom_fac));

        return ret; 
      }


    // general template that likes to blow up in your face
    template<typename A, typename B, typename C>
      data_store
      generate_tag(const rHandle<AbsRedstarInput_t> & snk, 
          const rHandle<AbsRedstarInput_t> & ins, 
          const rHandle<AbsRedstarInput_t> & src,
          const std::string &elem_id)
      {
        __builtin_trap(); 
      }    

    struct GENERATE_TAG_ERROR {};

    template<>
      data_store
      generate_tag<GENERATE_TAG_ERROR,GENERATE_TAG_ERROR,GENERATE_TAG_ERROR>
      (const rHandle<AbsRedstarInput_t> & snk, 
          const rHandle<AbsRedstarInput_t> & ins, 
          const rHandle<AbsRedstarInput_t> & src,
          const std::string &elem_id)
      {
        data_store dummy; 
        std::cerr << __PRETTY_FUNCTION__ << ": Error in this context" << std::endl; 
        __builtin_trap(); 
        return dummy; 
      } 


    // specialization 
    template<>
      data_store
      generate_tag<
      RedstarSingleParticleMesonInput,
      RedstarUnimprovedVectorCurrentInput,
      RedstarSingleParticleMesonInput>(const rHandle<AbsRedstarInput_t> & a, 
          const rHandle<AbsRedstarInput_t> & b, 
          const rHandle<AbsRedstarInput_t> & c,
          const std::string &elem_id)
      {
        data_store ret; 

        rHandle<RedstarSingleParticleMesonInput> snk(a);
        rHandle<RedstarUnimprovedVectorCurrentInput> ins(b);
        rHandle<RedstarSingleParticleMesonInput> src(c); 

        std::stringstream ss, mat_elem; 
        ss << snk->sname() << ".V_" << ins->lorentz 
          << "." << src->sname() << std::endl; 

        mat_elem << "_" << snk->H << "_" << src->H; 

        ret.file_id = ss.str(); 
        ret.jmu = ins->lorentz;
        ret.p_f = snk->mom;
        ret.p_i = src->mom; 
        ret.mat_elem_id = mat_elem.str(); 

        return ret;
      }

    // the next series of templates auto unrolls into the complicated 
    //     if else structure that is systemic across hadspec codes
    //     NB: unspecialized templates blow up!!

    // return value
    typedef data_store 
      (*FunctionPtr)(const rHandle<AbsRedstarInput_t> & ,
          const rHandle<AbsRedstarInput_t> & ,
          const rHandle<AbsRedstarInput_t> & ,
          const std::string &);

    // last guy in call chain 
    template<typename A, typename B, typename C> 
      FunctionPtr
      function_factory(void)
      {
        return &generate_tag<A,B,C>; 
      }

    // deal with the source
    template<typename A, typename B> 
      FunctionPtr
      function_factory(const rHandle<AbsRedstarInput_t> &  source)
      {
        FunctionPtr foo = &generate_tag<GENERATE_TAG_ERROR,
        GENERATE_TAG_ERROR,
        GENERATE_TAG_ERROR>;  // a blow up condition 

        if ( source->type() == Stringify<RedstarSingleParticleMesonInput>() ) 
        {
          foo = function_factory<A,B,RedstarSingleParticleMesonInput>(); 
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, source type " 
            << source->type() << " is unrecognized" << std::endl; 
        }

        return foo; 
      }

    // deals with the insertion 
    template<typename A> 
      FunctionPtr
      function_factory(const rHandle<AbsRedstarInput_t> & insertion, 
          const rHandle<AbsRedstarInput_t> & source)
      {
        FunctionPtr foo = &generate_tag<GENERATE_TAG_ERROR,
        GENERATE_TAG_ERROR,
        GENERATE_TAG_ERROR>;  // a blow up condition 

        if ( insertion->type() == Stringify<RedstarUnimprovedVectorCurrentInput>() ) 
        {
          foo = function_factory<A,RedstarUnimprovedVectorCurrentInput>(source); 
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, insertion type " 
            << insertion->type() << " is unrecognized" << std::endl; 
        }
        return foo; 
      }


    // deals with the sink 
    FunctionPtr 
      function_factory(const rHandle<AbsRedstarInput_t> & sink,
          const rHandle<AbsRedstarInput_t> & insertion,
          const rHandle<AbsRedstarInput_t> & source)
      {
        FunctionPtr foo = &generate_tag<GENERATE_TAG_ERROR,
        GENERATE_TAG_ERROR,
        GENERATE_TAG_ERROR>;  // a blow up condition 

        if ( sink->type() == Stringify<RedstarSingleParticleMesonInput>() ) 
        {
          foo = function_factory<RedstarSingleParticleMesonInput>(insertion,source); 
        }
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, sink type " 
            << sink->type() << " is unrecognized" << std::endl; 
        }

        return foo; 
      }


    // now put the pieces together, 
    //    1) check sanity
    //    2) get the data_store for this set of inputs 
    //           -- unspecialized templates are rigged 
    //              to blow up with gcc 
    //    3) make and return the tag
    LatticeMultiDataTag
      do_work(const std::vector<rHandle<AbsRedstarInput_t> > &elem, 
          const double mom_factor, 
          const double m_snk, 
          const double m_src, 
          const std::string &elem_id)
      {
        POW2_ASSERT(elem.size() == 3);
        FunctionPtr foo = function_factory(elem[0],elem[1],elem[2]); // this can blow 
        return generate_tag(
            (*foo)(elem[0],elem[1],elem[2],elem_id) // some data store -- can blow
            ,mom_factor,m_snk,m_src);  // the rest of the junk 
      }

  } // anonymous


  // the hard part
  LatticeMultiDataTag
    retrieve_tag(const std::vector<rHandle<AbsRedstarInput_t> > &elem, 
        const double mom_factor, 
        const double m_snk, 
        const double m_src, 
        const std::string &elem_id)
    {
      return do_work(elem,mom_factor,m_snk,m_src,elem_id); 
    }


  std::vector<TaggedEnsemRedstarNPtBlock>
    tag_lattice_xml(const AbstractMergeNamedObject * const ptr, 
        const double mom_factor, 
        const double m_snk, 
        const double m_src, 
        const std::string &elem_id)
    { 
      const AbsRedstarMergeNPtData_t * data_ptr = &(ptr->param->my_data); 
      unsigned int sz = data_ptr->npoint.size(); 
      POW2_ASSERT(sz == data_ptr->input.size()); 
      std::vector<TaggedEnsemRedstarNPtBlock> ret(sz);

      int idx,stop = int(sz); 
      
#ifdef PARALLEL_TAG_REDSTAR_DATA
#pragma omp parallel for shared(idx,stop) schedule(dynamic,100)
#endif 

      for (idx = 0; idx < stop; ++idx)
      {
        ret[idx] = TaggedEnsemRedstarNPtBlock(data_ptr->npoint[idx],
            retrieve_tag(data_ptr->input[idx],mom_factor,m_snk,m_src,elem_id)); 
      }
      
#ifdef PARALLEL_TAG_REDSTAR_DATA
#pragma omp barrier
#endif 
      return ret;
    }



} // radmat


