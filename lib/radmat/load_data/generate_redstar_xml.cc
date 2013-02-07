/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : generate_redstar_xml.cc

 * Purpose :

 * Creation Date : 03-12-2012

 * Last Modified : Thu Jan 31 15:57:39 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "generate_redstar_xml.h"
#include "radmat/load_data/invert_subduction.h"
#include "radmat/load_data/simple_world.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/pow2assert.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "hadron/hadron_npart_irrep.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "hadron/ensem_filenames.h"
#include "radmat_overlap_key_val_db.h"
#include "ensem/ensem.h"
#include "formfac/formfac_qsq.h"
#include "semble/semble_file_management.h"
#include "itpp/itbase.h"
#include <iostream>
#include <sstream>

#if 0



// struct to hold the xml and coefficients associated with a given lattice xml mat elem
struct redstarCircularMatElem_t
{
  typedef Hadron::KeyHadronNPartNPtCorr_t::NPoint_t redKey; 
  typedef RadmatExtendedKeyHadronNPartIrrep_t normKey;

  struct OperatorKeyData
  {
    OperatorKeyData(const redKey &red, const normKey &norm)
      : operatorKey(red) , normalizationKey(norm) {}

    redKey operatorKey;
    normKey normalizationKey; 
  };


  typedef ObjExpr_t<ENSEM::Complex,OperatorKeyData> subducedOp; 
  typedef ListObjExpr_t<ENSEM::Complex,subducedOp> listSubducedOp; 

  typedef ObjExpr_t<ENSEM::Complex,redKey> subducedInsertion;
  typedef ListObjExpr_t<ENSEM::Complex,subducedInsertion> listSubducedInsertion;


  // overload this constructor to deal with optimized operators???
  redstarCircularMatElem_t(const simpleWorld::ContinuumMatElem &, const std::string &source_id, const std::string &sink_id); 

  // the extra name of the state for the normalization database
  const std::string m_source_id;
  const std::string m_sink_id; 

  // lists of subduced operators and associated coefficients 
  listSubducedOp m_source;
  listSubducedOp m_sink
    std::pair<bool,listSubducedInsertion> m_time; 
  std::pair<bool,listSubducedInsertion> m_plus;
  std::pair<bool,listSubducedInsertion> m_zero;
  std::pair<bool,listSubducedInsertion> m_minus; 
};
#endif


// CIRCULAR BASIS STUFF
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

namespace radmat
{

  namespace 
  {
    // save some typing and in the process make this virtually unreadable to anyone other than myself..
    typedef redstarCircularMatElem_t::redKey redKey;
    typedef redstarCircularMatElem_t::normKey normKey;
    typedef redstarCircularMatElem_t::subducedOp subducedOp;
    typedef redstarCircularMatElem_t::listSubducedOp listSubducedOp;
    typedef redstarCircularMatElem_t::subducedInsertion subducedInsertion;
    typedef redstarCircularMatElem_t::listSubducedInsertion listSubducedInsertion;


    listSubducedOp getLatticeSubducedOp(const simpleWorld::ContinuumStatePrimitive &cont,
        const std::string &pid,
        const int t_slice)
    {
      listSubducedOp ret; 
      ContinuumBosonExprPrimitive boson(cont.J,cont.parity,cont.H,Hadron::generateLittleGroup(cont.mom)); 
      ListLatticeIrrepExpr_t lattice_boson = invertSubduction(boson); 
      ListLatticeIrrepExpr_t::const_iterator it;

      for(it = lattice_boson.begin(); it != lattice_boson.end(); ++it)
      {
        std::stringstream concat_name; 
        concat_name << cont.name << "_" << it->m_obj.irrep;

#if 0
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << cont.name << std::endl;

        if(boson.group == "Oh")
        {
          concat_name << "_" << it->m_obj.irrep;
        }
        else
        {
          std::cout << __PRETTY_FUNCTION__ << std::endl;
          std::cout << "group " << it->m_obj.group << std::endl;
          std::cout << "irrep " << it->m_obj.irrep << std::endl;

          concat_name << "_" << it->m_obj.irrep;
        }

        std::cout << __PRETTY_FUNCTION__ << std::endl;
        std::cout << concat_name.str() << std::endl; 
#endif

        Hadron::KeyHadronNPartIrrep_t base; 
        base.ops.resize(1); 
        base.ops[1].name = concat_name.str(); 
        base.ops[1].mom_type = FF::canonicalOrder(cont.mom); 
        base.row = it->m_obj.row; 
        base.twoI_z = cont.twoI_z; 
        base.mom = cont.mom;
        base.creation_op = cont.creation_op;
        base.smearedP = cont.smearedP;

        Hadron::KeyHadronNPartNPtCorr_t::NPoint_t npt;
        npt.t_slice = t_slice;
        npt.irrep = base;   

        RadmatExtendedKeyHadronNPartIrrep_t norm(pid,base);     
        redstarCircularMatElem_t::OperatorKeyData key(npt,norm);

        ret.push_back(subducedOp(it->m_coeff,key));
      }   
      return ret; 
    }


    listSubducedOp getLatticeSubducedOp(const simpleWorld::ContinuumMatElem::State &state, const std::string &pid)
    {
      return getLatticeSubducedOp(state.state,pid,state.t_slice);
    }


    struct circLorentzInsertion
    {
      circLorentzInsertion(void)
        : time(std::pair<bool,listSubducedInsertion>(false,listSubducedInsertion())) , plus(time), zero(time) , minus(time) {}

      std::pair<bool,listSubducedInsertion> time;
      std::pair<bool,listSubducedInsertion> plus;
      std::pair<bool,listSubducedInsertion> zero;
      std::pair<bool,listSubducedInsertion> minus; 
    };


    circLorentzInsertion getLatticeSubducedInsertion(const simpleWorld::ContinuumInsertion &ins, const std::string &pid)
    {
      circLorentzInsertion ret;
      listSubducedOp work;
      listSubducedOp::const_iterator work_it;
      std::map<std::string,simpleWorld::ContinuumInsertion::op_insertion>::const_iterator it;
      ADATXML::Array<int> q;
      bool haveq(false); 

      for(it = ins.insertion_map.begin(); it != ins.insertion_map.end(); ++it)
      {
        work = getLatticeSubducedOp(it->second,pid,ins.t_slice);

        if(!!!haveq)
        {
          q = it->second.mom; 
          haveq = true; 
        }

        listSubducedInsertion ins_list;

        for(work_it = work.begin(); work_it != work.end(); ++work_it)
          ins_list.push_back(subducedInsertion(work_it->m_coeff,work_it->m_obj.operatorKey));

        if(it->first == "t")
          ret.time = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else if(it->first == "p")
          ret.plus = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else if(it->first == "0")
          ret.zero = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else if(it->first == "m")
          ret.minus = std::pair<bool,listSubducedInsertion>(true,ins_list);
        else
        {
          std::cout << __func__ << ": error: can't match to a circular lorentz index" << std::endl;
          exit(1);
        }

      }

      return ret; 
    }


    std::string get_string_list(const redstarCartMatElem::redstarCartMatElemLorentzComponent &comp)
    {
      redstarCartMatElem::redstarCartMatElemLorentzComponent::const_iterator it;
      std::stringstream ss; 

      std::map<std::string,ENSEM::Complex> collapse_xml;
      std::map<std::string,ENSEM::Complex>::iterator map_it; 

      for(it = comp.begin(); it != comp.end(); ++it)
      {
        std::string name =  Hadron::ensemFileName(it->m_obj.redstar_xml) ;
        map_it = collapse_xml.find(name);

        if(map_it != collapse_xml.end())
        {
          ENSEM::Complex tmp = it->m_coeff + map_it->second; 
          map_it->second = tmp;
        }
        else
        {
          collapse_xml.insert(std::pair<std::string,ENSEM::Complex>(name,it->m_coeff));
        }
      }


      for(map_it = collapse_xml.begin(); map_it != collapse_xml.end(); ++map_it)
        ss << SEMBLE::toScalar(map_it->second) << " * " << map_it->first << std::endl;

      return ss.str();
    }




    //////////////////////////
    //////////////////////////
    //////////////////////////
    //////////////////////////


    // a bunch of overloads and functions to move to cartesian coordinates


    const std::complex<double> one(std::complex<double>(1.,0.));
    const std::complex<double> cplx_i(std::complex<double>(0.,1.));
    const double root2(sqrt(2.));
    const double root2inv(1./sqrt(2.)); 


    redstarCircularMatElem_t::listSubducedInsertion 
      operator*(const std::complex<double> &c, const redstarCircularMatElem_t::listSubducedInsertion &l)
      {
        return SEMBLE::toScalar(c)*l; 
      }


    itpp::Vec<redstarCircularMatElem_t::listSubducedInsertion> 
      operator*(const itpp::Mat<std::complex<double> > &m, const itpp::Vec<redstarCircularMatElem_t::listSubducedInsertion> &v)
      {

        POW2_ASSERT(v.size() == m.cols()); 

        itpp::Vec<redstarCircularMatElem_t::listSubducedInsertion> ret(m.rows()); 
        redstarCircularMatElem_t::listSubducedInsertion tmp; 

        for(int row = 0; row < m.rows(); ++row)
          for(int col = 0; col < v.size(); ++col)
          {
            tmp = ret[row] + m(row,col)*v(col);
            ret[row] = tmp; 
          }

        return ret; 
      }

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
            if(itpp::round_to_zero(m(row,col), 0.00001) == std::complex<double>(0.,0.))
              continue; 

            ret(row) &= v(col);
          }

        return ret; 
      }

    itpp::Vec<std::complex<double> > eps3_z(const std::string &qn)
    {
      itpp::Vec<std::complex<double> > ret(3);
      ret.zeros(); 

      if(qn == "p")
      {
        ret[0] = -root2inv;
        ret[1] = -cplx_i/root2; 
      }
      else if(qn == "m")
      {
        ret[0] = root2inv;
        ret[1] = -cplx_i/root2;
      }
      else if(qn == "0")
      {
        ret[2] = 1;
      }
      else
      {
        std::cerr << "qn(" << qn << ") is ill defined, use p,m,0" << std::endl;
        exit(1);
      }

      return ret;
    }


    // columns are indexed -m to m
    // rows go -lambda to lambda

    itpp::Mat<std::complex<double> > Wigner_D(const ADATXML::Array<int> mom, const int J)
    {

      if((mom[0] == 0) &&(mom[1] == 0) && (mom[2] == 0))
      {
        itpp::Vec<std::complex<double> > one(3);
        one.ones();
        return itpp::diag(one);
      }



      Hadron::CubicCanonicalRotation_t rot = Hadron::cubicCanonicalRotation(mom);
      itpp::Mat<std::complex<double> > ret(2*J + 1, 2*J+1); 
      std::complex<double> tmp; 
      std::complex<double> zero(0.,0.); 

      const int tJ = 2*J;

      for(int lambda = -J; lambda < J +1; ++lambda)
        for(int m = -J; m < J +1; ++m)
        {
          tmp = SEMBLE::toScalar(Hadron::Wigner_D(tJ,2*m,2*lambda,rot.alpha,rot.beta,rot.gamma));

          if( itpp::round_to_zero(tmp,0.00001) == zero )
            ret(lambda+J,m+J) = zero;
          else
            ret(lambda+J,m+J) = tmp; 
        }

      return ret; 
    }



    // a bunch of hardwires since itpp and ensem don't get along..
    itpp::Mat<std::complex<double> > inver2Cart(const ADATXML::Array<int> mom)
    { 
      itpp::Mat<std::complex<double> > epsz(3,3), D;
      epsz.set_row(0,eps3_z("m"));
      epsz.set_row(1,eps3_z("0"));
      epsz.set_row(2,eps3_z("p"));

      D = Wigner_D(mom,1); 
      return itpp::inv(D*epsz); 
    }





    //////////////////////////
    //////////////////////////
    //////////////////////////
    //////////////////////////



    // the lattice matrix elements are in a helicity basis but the linear system
    // solvers are set up to use a cartesian basis.. use this to do a trivial change of basis
    struct redstarCartesianMatElem_p
    {
      // save some typing and in the process make this virtually unreadable to anyone other than myself..
      typedef redstarCircularMatElem_t::redKey redKey;
      typedef redstarCircularMatElem_t::normKey normKey;
      typedef redstarCircularMatElem_t::subducedOp subducedOp;
      typedef redstarCircularMatElem_t::listSubducedOp listSubducedOp;
      typedef redstarCircularMatElem_t::subducedInsertion subducedInsertion;
      typedef redstarCircularMatElem_t::listSubducedInsertion listSubducedInsertion;


      redstarCartesianMatElem_p(const simpleWorld::ContinuumMatElem &a, 
          const std::string &source_id, 
          const std::string &sink_id)
      {
        redstarCircularMatElem_t dum(a,source_id,sink_id);
        ensemble = a.ensemble; 
        m_source = dum.m_source;
        m_sink = dum.m_sink;
        makeCartesian(dum);
      }

      // move from a circular basis to the cartesian basis 
      void makeCartesian(const redstarCircularMatElem_t &tmp);

      ADATXML::Array<int> mom(void)
      {
        ADATXML::Array<int> pf,pi,q;
        pf = m_sink.begin()->m_obj.operatorKey.irrep.mom;
        pi = m_source.begin()->m_obj.operatorKey.irrep.mom;

        q.resize(3);
        q[0] = pi[0] - pf[0];
        q[1] = pi[1] - pf[1];
        q[2] = pi[2] - pf[2];

        return q; 
      }


      listSubducedOp m_source;
      listSubducedOp m_sink;
      std::pair<bool,listSubducedInsertion> m_t; 
      std::pair<bool,listSubducedInsertion> m_x;
      std::pair<bool,listSubducedInsertion> m_y;
      std::pair<bool,listSubducedInsertion> m_z;
      std::string ensemble;  
    };




    void redstarCartesianMatElem_p::makeCartesian(const redstarCircularMatElem_t &tmp)
    {
      itpp::Vec<listSubducedInsertion> j_lambda(3),j_i(3);
      itpp::Vec<bool> j_lambda_have(3),j_i_have(3);

      j_lambda[0] = tmp.m_minus.second;
      j_lambda[1] = tmp.m_zero.second; 
      j_lambda[2] = tmp.m_plus.second; 

      j_lambda_have[0] = tmp.m_minus.first;
      j_lambda_have[1] = tmp.m_zero.first; 
      j_lambda_have[2] = tmp.m_plus.first; 


      itpp::Mat<std::complex<double> > tform = inver2Cart(mom());

      j_i = tform * j_lambda;
      j_i_have = tform * j_lambda_have; 

      m_t = tmp.m_time;
      m_x = std::pair<bool,listSubducedInsertion>(j_i_have[0],j_i[0]);    
      m_y = std::pair<bool,listSubducedInsertion>(j_i_have[1],j_i[1]);  
      m_z = std::pair<bool,listSubducedInsertion>(j_i_have[2],j_i[2]);  

      /*      
              std::cout << "false = " << false << std::endl;
              std::cout << "j_lambda_have 0 , 1 , 2 ---" << j_lambda_have[0] 
              << "  " << j_lambda_have[1]
              << "  " << j_lambda_have[2] << std::endl;
              std::cout << "j_i_have 0 , 1 , 2 --- " << j_i_have[0]
              << "  " << j_i_have[1] 
              << "  " << j_i_have[2] << std::endl; 
       */
    }


  } // namespace anonomyous 



  redstarCircularMatElem_t::redstarCircularMatElem_t(const simpleWorld::ContinuumMatElem & cont,
      const std::string &source_id, 
      const std::string &sink_id)
    : m_source_id(source_id), m_sink_id(sink_id)
  {
    m_source = getLatticeSubducedOp(cont.source,m_source_id);
    m_sink = getLatticeSubducedOp(cont.sink,m_sink_id);

    /*
       std::cout << __func__ << std::endl;
       std::cout << "source " << cont.source << std::endl;
       std::cout << "sink   " << cont.sink << std::endl;
     */

    circLorentzInsertion foo = getLatticeSubducedInsertion(cont.insertion,std::string("foobar"));
    m_time = foo.time;
    m_plus = foo.plus;
    m_zero = foo.zero;
    m_minus = foo.minus;
  }


} // namespace radmat




// CARTESIAN BASIS STUFF
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

namespace radmat
{

  namespace
  {
    redstarCartMatElem::redstarCartMatElemLorentzComponent::threePointKey genThreePointKey(const listSubducedOp::ListObj_t &source, 
        const listSubducedInsertion::ListObj_t &ins,
        const listSubducedOp::ListObj_t &sink,
        const std::string &ensemble)
    {
      Hadron::KeyHadronNPartNPtCorr_t redstar_xml;
      redstar_xml.npoint.resize(3); // 1 based array
      redstar_xml.npoint[1] = sink.m_obj.operatorKey;
      redstar_xml.npoint[2] = ins.m_obj;
      redstar_xml.npoint[3] = source.m_obj.operatorKey;
      redstar_xml.ensemble = ensemble; 

      /*
         std::cout << __func__ << Hadron::ensemFileName(redstar_xml) << std::endl;
         std::cout << "        source " << source.m_obj.normalizationKey << std::endl;
         std::cout << "        sink   " << sink.m_obj.normalizationKey << std::endl;
       */

      return redstarCartMatElem::redstarCartMatElemLorentzComponent::threePointKey(redstar_xml,source.m_obj.normalizationKey,sink.m_obj.normalizationKey); 
    }

  } // namespace anonomyous


  redstarCartMatElem::redstarCartMatElemLorentzComponent::redstarCartMatElemLorentzComponent(const listSubducedOp &source,
      const std::pair<bool,listSubducedInsertion> &ins,
      const listSubducedOp &sink,
      const std::string &ensemble) : active(false)
  {

    // leave the list empty if the init variable from the xml was false for whatever reason
    if(ins.first)
    {
      listSubducedOp::const_iterator it_source,it_sink;
      listSubducedInsertion::const_iterator it_ins;

      // run the triple for loop and keep track of the subduction coefficients 
      for(it_source = source.begin(); it_source != source.end(); ++it_source)
        for(it_ins = ins.second.begin(); it_ins != ins.second.end(); ++it_ins)
          for(it_sink = sink.begin(); it_sink != sink.end(); ++it_sink)
            m_list.push_back(threePointXMLKey(it_source->m_coeff * it_ins->m_coeff * it_sink->m_coeff,
                  genThreePointKey(*it_source,*it_ins,*it_sink,ensemble)
                  )
                );
      active = true; 
    }

  }


  redstarCartMatElem::redstarCartMatElem(const simpleWorld::ContinuumMatElem &cont,
      const std::string &source_id,
      const std::string &sink_id)
    : m_elem(cont)
  {
    redstarCartesianMatElem_p foobar(cont,source_id,sink_id);

    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(0,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_t,foobar.m_sink,foobar.ensemble))
        );
    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(1,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_x,foobar.m_sink,foobar.ensemble))
        );
    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(2,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_y,foobar.m_sink,foobar.ensemble))
        );
    lorentz_components.insert(
        std::map<int,redstarCartMatElemLorentzComponent>::value_type(3,
          redstarCartMatElemLorentzComponent(foobar.m_source,foobar.m_z,foobar.m_sink,foobar.ensemble))
        );

    if(foobar.m_t.first)
      lorentz_components[0].active = true;
    if(foobar.m_x.first)
      lorentz_components[1].active = true;
    if(foobar.m_y.first)
      lorentz_components[2].active = true;
    if(foobar.m_z.first)
      lorentz_components[3].active = true;

  }


  redstarCartMatElem::~redstarCartMatElem(void)
  {
    std::map<int,redstarCartMatElemLorentzComponent>::const_iterator it; 

    for(it = lorentz_components.begin(); it != lorentz_components.end(); ++it)
    {
      std::string eqn = get_string_list(it->second);
      int idx = it->first; 

      std::string source,sink;
      source = simpleWorld::stateFileName(m_elem.source);
      sink = simpleWorld::stateFileName(m_elem.sink);
      std::stringstream ss;

      ss << sink << "___ins_" << idx <<"___" << source;

      std::string pth = SEMBLE::SEMBLEIO::getPath();
      std::stringstream output;
      output << pth << "cont_projected_matrix_elements";
      SEMBLE::SEMBLEIO::makeDirectoryPath(output.str());
      output << "/xml";
      SEMBLE::SEMBLEIO::makeDirectoryPath(output.str());
      output << "/" << ss.str() << ".redstar.xml";

      std::ofstream out(output.str().c_str());
      out << eqn << std::endl;
      out.close();

    }



  }



} // namespace radmat

