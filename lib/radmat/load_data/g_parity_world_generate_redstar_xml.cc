/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : g_parity_world_generate_redstar_xml.cc

 * Purpose :

 * Creation Date : 25-04-2013

 * Last Modified : Thu May  9 16:23:42 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "g_parity_world_generate_redstar_xml.h"
#include "invert_subduction.h"

#include "semble/semble_meta.h"

#include "adat/map_obj.h"
#include "ensem/ensem.h"
#include "hadron/ensem_filenames.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "formfac/formfac_qsq.h"

#include "radmat/utils/polarisation_tensors.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/obj_expr_t.h"

#include "itpp/itbase.h"

#include <map>

namespace radmat
{


  //////////////////////////
  //////////////////////////
  //////////////////////////
  //////////////////////////


  // a bunch of overloads and functions to move to cartesian coordinates

  namespace
  {


    typedef generateCircularRedstarXML::listNPointKey listNPointKey;


    // a debugging method
    void screen_dump (const listNPointKey &k)
    {

#if 0
      std::cout << __func__ << ": nele = " << k.m_expr.size() << std::endl;

      //     for(int i = 0; i < k.m_expr.size(); ++i)
      //       std::cout << ensemFileName( k.m_expr[i].m_obj ) << std::endl;

      listNPointKey::const_iterator it;
      for(it = k.begin(); it != k.end(); ++it)
      {
        std::cout << it->m_obj << std::endl;
        std::cout << ensemFileName(it->m_obj) << std::endl;
      }

#endif
    } 


    listNPointKey operator*(const std::complex<double> &c, const listNPointKey &l)
    {
      return SEMBLE::toScalar(c)*l; 
    }

    itpp::Vec<listNPointKey> operator*(const itpp::Mat<std::complex<double> > &m,
        const itpp::Vec<listNPointKey> &v)
    {
      POW2_ASSERT(v.size() == m.cols()); 
      itpp::Vec<listNPointKey> ret(m.rows()); 

      for(int row = 0; row < m.rows(); ++row)
        for(int col = 0; col < m.cols(); ++col)
        {
          if(std::norm(m(row,col)) > 0.0001)
          {
            // init 
            if(ret[row].m_expr.size() == 0)
              ret[row] =  m(row,col)*v(col);
            else
              ret[row] = ret[row] + m(row,col)*v(col);
          }
        }

      return ret; 
    }


    // multiply a matrix against a vector of bools
    itpp::Vec<bool>
      operator*(const itpp::Mat<std::complex<double> > &m, const itpp::Vec<bool> &v)
      {
        POW2_ASSERT(v.size() == m.cols()); 
        itpp::Vec<bool> ret(m.rows());

        for(int row = 0; row < ret.size(); ++row)
          ret(row) = true; 


        for(int row = 0; row < m.rows(); ++row)
          for(int col = 0; col < v.size(); ++col)
          {
            // phases may not cancel exactly but in this context 0.000001 is zero
            if(std::norm(m(row,col)) > 0.0001)
              continue; 

            ret(row) &= v(col);
          }

        return ret; 
      }


    // rows correspond to a helicity (+ 0 -)
    // cols are the cartesian coordinate (x,y,z)
    // then M(row,col) = epsilon_cart(lambda)  -->  j^lambda = M * j^cartesian 
    itpp::Mat<std::complex<double> > eps3d(const ADATXML::Array<int> &mom , const bool create)
    {
      Tensor<std::complex<double>, 1 > tmp;
      genPolTens3D<1> eps(mom);
      itpp::Mat<std::complex<double> > eps3(3,3); 

      for(int h = 1; h > -2; --h)
      {
        tmp = eps.get(h);

        for(int i = 0; i < 3; ++i)
          if(create)
            eps3(1-h,i) = tmp[i];
          else
            eps3(1-h,i) = std::conj(tmp[i]); 

      }

      return eps3; 
    }


    // invert the matrix eps
    itpp::Mat<std::complex<double> > invert2Cart(const ADATXML::Array<int> mom, const bool create)
    { 
      return itpp::round_to_zero(itpp::inv(eps3d(mom,create)),0.00001);
    }

    std::string stringy_mom(const ADATXML::Array<int> mom)
    {
      ADATXML::Array<int> can = FF::canonicalOrder(mom);
      std::stringstream ss;
      ss << "p" << can[0] << can[1] << can[2];
      return ss.str(); 
    }

  } // anonomyous 




  void redstarSubducedOperator::generate_subduction_list(const gParityWorld::GParityState& e,
      const int t_slice)
  {
    orig = e; 
    ContinuumBosonExprPrimitive meson(e.J,e.parity,e.H,Hadron::generateLittleGroup(e.mom)); 

    ListLatticeIrrepExpr_t lattice_meson;

    if(e.creation_op)
      lattice_meson = invertSubduction(meson); 
    else
      lattice_meson = conj(invertSubduction(meson)); 

    ListLatticeIrrepExpr_t::const_iterator it; 
    for(it = lattice_meson.begin(); it != lattice_meson.end(); ++it)
    {
      std::stringstream name;
      name << e.name << "_";
      if(e.isProjected)
        name << stringy_mom(e.mom) << "_";
      name << it->m_obj.irrep; 

      Hadron::KeyHadronNPartIrrep_t base; 
      base.ops.resize(1); 
      base.ops[1].name = name.str(); 
      base.ops[1].mom_type = FF::canonicalOrder(e.mom); 
      base.row = it->m_obj.row;
      base.twoI_z = e.twoI_z;
      base.mom = e.mom;
      base.creation_op = e.creation_op;
      base.smearedP = e.smearedP;

      redKey npt; 
      npt.t_slice = t_slice; 
      npt.irrep = base; 

      subduced = subduced + irrepOperator(it->m_coeff,npt); 
    }
  }

  // -- 

  void redstarSubducedState::generate_subduction_list(
      const gParityWorld::GParityHelicityMatrixElement::State &e)
  {
    redstarSubducedOperator work; 
    work.generate_subduction_list(e.state,e.t_slice); 
    subduced = work.subduced;
    orig = e; 
  }

  // -- 




  void redstarSubducedPhoton::generate_subduction_list(
      const gParityWorld::GParityInsertion::photon &p, const int t_slice)
  {
    gParityWorld::GParityInsertion::photon::const_iterator it; 

    for(it = p.begin(); it != p.end(); ++it)
    {
      redstarSubducedOperator work;
      work.generate_subduction_list(it->m_obj,t_slice); 

      // the charge coefficient sitting in front of the quark is real but the 
      // subduction coefficients are complex, promote to 'complex' charge here
      subduced = subduced + SEMBLE::toScalar(std::complex<double>(it->m_coeff,0.)) * work.subduced;
    }
    orig = p; 
  }

  // -- 

  void redstarSubducedHelicityInsertion::generate_subduction_list(
      const gParityWorld::GParityInsertion &e)
  {
    gParityWorld::GParityInsertion::map_t::const_iterator it; 

    for(it = e.begin(); it != e.end(); ++it)
    {
      redstarSubducedPhoton work; 
      work.generate_subduction_list(it->second,e.t_slice);
      insertion_map.insert(map_t::value_type(it->first,work.subduced)); 
    }
    orig = e; 
  }

  // -- 

  Hadron::KeyHadronNPartNPtCorr_t 
    makeNPartNPoint(const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &sink,
        const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &insertion,
        const Hadron::KeyHadronNPartNPtCorr_t::NPoint_t &source,
        const std::string &ensemble)
    {
      Hadron::KeyHadronNPartNPtCorr_t ret;
      ret.npoint.resize(3); 
      ret.npoint[1] = sink;
      ret.npoint[2] = insertion;
      ret.npoint[3] = source;
      ret.ensemble = ensemble; 
      return ret; 
    } 

  // -- 

  void mergeSubducedOperators3pt::merge_subduction_lists(const listIrrepOperator &sinks,
      const listIrrepOperator &insertions,
      const listIrrepOperator &sources,
      const std::string &ensemble)
  {
    m_sink = sinks;
    m_insertion = insertions;
    m_source = sources; 

    std::map<std::string,std::pair<ENSEM::Complex,NPointKey> > coeff_map; 

    listIrrepOperator::const_iterator sink,insertion,source;

    for(sink = sinks.begin(); sink != sinks.end(); ++sink)
      for(insertion = insertions.begin(); insertion != insertions.end(); ++insertion)
        for(source = sources.begin(); source != sources.end(); ++source)
        {
          ENSEM::Complex coeff = sink->m_coeff*insertion->m_coeff*source->m_coeff;
          NPointKey npt = makeNPartNPoint(sink->m_obj,insertion->m_obj,source->m_obj,ensemble); 
          std::string key = ensemFileName(npt); 

          if(coeff_map.find(key) != coeff_map.end())
            coeff += coeff_map[key].first;

          // reinsert
          coeff_map[key] = std::pair<ENSEM::Complex,NPointKey>(coeff,npt); 
        }

    // now make the list and include the guys with zero coefficients 
    // which we will deal with at a later pass
    std::map<std::string,std::pair<ENSEM::Complex,NPointKey> >::const_iterator it; 

    for(it = coeff_map.begin(); it != coeff_map.end(); ++it)
      subduced = subduced + objNPointKey(it->second.first,it->second.second); 
  }

  // -- 

  void generateCircularRedstarXML::fill(const std::string &id, const map_t& m, data_t &d, 
      const listIrrepOperator &sink, const listIrrepOperator &source,
      const std::string &ensemble)
  {
    mergeSubducedOperators3pt work;
    if(m.find(id) != m.end())
    {
      work.merge_subduction_lists(sink,m.find(id)->second,source,ensemble);
      d = data_t(true,work.subduced); 
    }
    else
    {
      //     std::cout << __func__ << ": missed id = " << id << std::endl;
      d = data_t(false,work.subduced); // empty list
    }
  }


  void generateCircularRedstarXML::generate_subduction_list(
      const gParityWorld::GParityHelicityMatrixElement &e)
  {
    redstarSubducedHelicityInsertion insertion;
    redstarSubducedState sink, source; 

    insertion.generate_subduction_list(e.insertion);
    sink.generate_subduction_list(e.sink);
    source.generate_subduction_list(e.source);

    // fill 'er up      
    fill("t", insertion.insertion_map,time,sink.subduced,source.subduced,e.ensemble); 
    fill("p", insertion.insertion_map,plus,sink.subduced,source.subduced,e.ensemble); 
    fill("0", insertion.insertion_map,zero,sink.subduced,source.subduced,e.ensemble); 
    fill("m", insertion.insertion_map,minus,sink.subduced,source.subduced,e.ensemble);

    orig = e; 
  }

  // --

  generateCartesianRedstarXML::data_t 
    generateCartesianRedstarXML::combine_duplicates(const generateCartesianRedstarXML::data_t &d) 
    {
      generateCartesianRedstarXML::data_t ret;
      ret.first = d.first; 

      // break on fail with an empty list
      if(!!! ret.first)
        return ret; 

      // loop and combine coefficients for duplicate npoints
      std::map<std::string,
        std::pair<ENSEM::Complex,generateCartesianRedstarXML::NPointKey> > coeff_map;
      generateCartesianRedstarXML::listNPointKey::const_iterator it;

      for(it = d.second.begin(); it != d.second.end(); ++it)
      {
        ENSEM::Complex coeff = it->m_coeff;

        std::string key = ensemFileName(it->m_obj); 

        if(coeff_map.find(key) != coeff_map.end())
          coeff += coeff_map[key].first;

        // reinsert
        coeff_map[key] = 
          std::pair<ENSEM::Complex,generateCartesianRedstarXML::NPointKey>(coeff,it->m_obj); 
      }

      // write the output list
      std::map<std::string,
        std::pair<ENSEM::Complex,generateCartesianRedstarXML::NPointKey> >::const_iterator map_it;
      for(map_it = coeff_map.begin(); map_it != coeff_map.end(); ++map_it)
        ret.second = ret.second 
          + generateCartesianRedstarXML::objNPointKey(map_it->second.first,map_it->second.second);

      return ret; 
    } 

  void generateCartesianRedstarXML::generate_subduction_list(
      const gParityWorld::GParityHelicityMatrixElement &e)
  {

    generateCircularRedstarXML tmp; 
    tmp.generate_subduction_list(e); 

    // time is simple
    t = tmp.time;

    // space is not hard  
    ADATXML::Array<int> q = e.insertion.mom;
    itpp::Mat<std::complex<double> > M = invert2Cart(q,e.insertion.create()); 

    itpp::Vec<bool> circb(3), cartb;
    circb[0] = tmp.plus.first;
    circb[1] = tmp.zero.first;
    circb[2] = tmp.minus.first;

    itpp::Vec<listNPointKey> circl(3), cartl;
    circl[0] = tmp.plus.second;
    circl[1] = tmp.zero.second;
    circl[2] = tmp.minus.second;

    cartb = M * circb; 
    cartl = M * circl; 

#if 0
    std::cout << __func__ << std::endl;

    std::cout << M << std::endl;

    //  std::cout << e << std::endl;

    std::cout << " pos " << tmp.plus.first << std::endl;
    screen_dump( circl[0] );

    std::cout << " zero " << tmp.zero.first << std::endl;
    screen_dump( circl[1] );

    std::cout << " neg " << tmp.minus.first << std::endl;
    screen_dump( circl[2] );
#endif

    x = combine_duplicates(data_t(cartb[0],cartl[0]));
    y = combine_duplicates(data_t(cartb[1],cartl[1]));
    z = combine_duplicates(data_t(cartb[2],cartl[2])); 

    orig = e; 
  }



}
