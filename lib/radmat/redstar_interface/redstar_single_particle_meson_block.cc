/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_single_particle_meson_block.cc

 * Purpose :

 * Creation Date : 11-11-2013

 * Last Modified : Thu 14 Nov 2013 03:09:48 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_single_particle_meson_block.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/construct_data/invert_subduction.h"
#include "radmat/utils/polarisation_tensors.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/obj_expr_t.h"

#include "adat/map_obj.h"
#include "ensem/ensem.h"
#include "hadron/ensem_filenames.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "formfac/formfac_qsq.h"
#include "itpp/itbase.h"

#include <exception>
#include <iostream>
#include <sstream>


#define  DEBUG_MSG_ON
#define  DEBUG_HANDLE_ON
#include "debug_props.h"

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


    // stringy mom
    std::string stringmom(const ADATXML::Array<int> &inp)
    {
      ADATXML::Array<int> can = FF::canonicalOrder(inp);
      std::stringstream ss;
      ss << can[0] << can[1] << can[2];
      return ss.str(); 
    }

    std::string stringy_mom(const ADATXML::Array<int> mom)
    {
      ADATXML::Array<int> can = FF::canonicalOrder(mom);
      std::stringstream ss;
      ss << "p" << can[0] << can[1] << can[2];
      return ss.str(); 
    }

    std::string string_mom_no_space(const ADATXML::Array<int> &b)
    {
      std::stringstream ss; 
      for(int i = 0; i < b.size(); ++i)
        ss << b[i];
      return ss.str(); 
    }

    // get the star
    ADATXML::Array<ADATXML::Array<int> > star_p(const ADATXML::Array<ADATXML::Array<int> > &pin)
    {

      std::list<ADATXML::Array<int> > work;
      std::map<std::string,std::list<ADATXML::Array<int> > > seen;
      std::map<std::string,std::list<ADATXML::Array<int> > >::const_iterator it;
      std::list<ADATXML::Array<int> >::const_iterator listit; 

      int nele(0); 

      const unsigned int sz = pin.size(); 

      for(unsigned int elem = 0; elem < sz; ++elem)
      {
        std::string p = stringmom(pin[elem]);
        it = seen.find(p);
        if (it == seen.end())
        {
          work = Hadron::generateLittleGroupMom(Hadron::generateLittleGroup(pin[elem]),
              FF::canonicalOrder(pin[elem]));
          nele += work.size();
          seen.insert(std::pair<std::string,std::list<ADATXML::Array<int> > >(p,work)); 
        }        
      }


      ADATXML::Array<ADATXML::Array<int> > ret(nele);
      nele = 0;

      for(it = seen.begin(); it != seen.end(); ++it)
        for(listit = it->second.begin(); listit != it->second.end(); ++listit)
        {
          ret[nele] = *listit; 
          ++nele; 
        }

      return ret; 
    }



    EnsemRedstarBlock handle_work(const RedstarSingleParticleMesonInput &e)
    {
      DEBUG_MSG(entering);
      // return value
      EnsemRedstarBlock ret; 

      // play the subduction dance
      ContinuumBosonExprPrimitive meson(e.J,e.parity,e.H,Hadron::generateLittleGroup(e.mom)); 

      // the lattice version 
      ListLatticeIrrepExpr_t lattice_meson;

      // handle conventions
      if(e.creation_op)
        lattice_meson = invertSubduction(meson); 
      else
        lattice_meson = conj(invertSubduction(meson)); 

      DEBUG_MSG(inverted subduction);

      // do this sum
      // O_{J,H} \sim \sum_{\lambda,\mu} S_{J,H}^{\lambda,\mu} O^{\lambda,\mu}
      ListLatticeIrrepExpr_t::const_iterator it; 
      for(it = lattice_meson.begin(); it != lattice_meson.end(); ++it)
      {
        std::stringstream name;
        name << e.name << "_";
        if(e.isProjected)
          name << stringy_mom(e.mom) << "_";
        name << it->m_obj.irrep; 

        Hadron::KeyHadronNPartIrrep_t base; 

        //    IN WITH THAT NEW STUFF
        base.op.ops.resize(1); 
        base.op.ops[1].name = name.str(); 
        base.op.ops[1].mom_type = FF::canonicalOrder(e.mom); 

        base.row = it->m_obj.row;
        base.twoI_z = e.twoI_z;
        base.mom = e.mom;
        base.creation_op = e.creation_op;
        base.smearedP = e.smearedP;

        Hadron::KeyHadronNPartNPtCorr_t::NPoint_t npt; 
        npt.t_slice = e.t_slice; 
        npt.irrep = base; 

        ret = ret + EnsemRedstarBlock::ListObj_t(it->m_coeff,npt); 
      }
      DEBUG_MSG(exiting);
      return ret; 
    }

    std::string doPrint(const ADATXML::Array<int> &t)
    {
      std::stringstream ss; 
      for(int i = 0; i < t.size(); ++i)
        ss << " " << t[i];
      return ss.str();  
    } 

    std::string doPrint(const ADATXML::Array< ADATXML::Array<int> > &t)
    {
      std::stringstream ss; 
      for(int i = 0; i < t.size(); ++i)
        ss << doPrint( t[i] ) << std::endl;
      return ss.str();  
    } 


  } // anonomyous 

  std::string 
    RedstarSingleParticleMesonInput::write(void) const
    {
      std::stringstream ss; 
      ss << "J= " << J << " H= " << H << " par= " << parity
        << " mom= " << doPrint(mom) << " twoI_z= " << twoI_z << " name = " << name 
        << " creation_op= " << creation_op << " smearedP= " << smearedP
        << " isProjected= " << isProjected << " t_slice= " << t_slice;
      return ss.str(); 
    }


  std::string 
    RedstarSingleParticleMesonInput::sname(void) const
    {
      std::stringstream ss; 
      ss << name << "_p" << string_mom_no_space(FF::canonicalOrder(mom))
        << ",J" << J << ",H" << H 
        << ",p" << string_mom_no_space(mom)
        << ",Iz" << twoI_z;
      return ss.str(); 
    }

  EnsemRedstarBlock
    RedstarSingleParticleMesonBlock::operator()(const AbsRedstarInput_t * base) const
    {
      RedstarSingleParticleMesonInput dummy; 
      POW2_ASSERT(base->type() == dummy.type());   

      try
      {
        const RedstarSingleParticleMesonInput *derived ; 
        derived = dynamic_cast<const RedstarSingleParticleMesonInput*>(base); 
        dummy = *derived; 
      }
      catch(std::exception &e)
      {
        std::cout << __PRETTY_FUNCTION__ <<  e.what() << std::endl; 
      }

      return handle_work(dummy); 
    }



  std::string 
    RedstarSingleParticleMesonXML::write(void) const
    {
      std::stringstream ss; 
      ss << "J= " << J << " H= " << doPrint(H) << " par= " << parity
        << " fill_star= " << fill_star 
        << " twoI_z= " << twoI_z << " name = " << name 
        << " creation_op= " << creation_op << " smearedP= " << smearedP
        << " isProjected= " << isProjected << " t_slice= " << t_slice
        << " momlist: \n " << doPrint(mom); 
      return ss.str(); 
    }

  void 
    RedstarSingleParticleMesonXML::read(ADATXML::XMLReader &xml, const std::string &path)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"J",J,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"H",H,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"parity",parity,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"fill_star",fill_star,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"mom",mom,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"twoI_z",twoI_z,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"name",name,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"creation_op",creation_op,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"smearedP",smearedP,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"isProjected",isProjected,__PRETTY_FUNCTION__);   
      doXMLRead(ptop,"t_slice",t_slice,__PRETTY_FUNCTION__);   


      if ( H.size() == 0) 
      {
        H.resize(2*J + 1);
        for(int h = -J; h < J+1; ++h)
          H[h+J] = h; 
      }

      if(fill_star)
        mom = star_p(mom); 


      objFunctorPtr = new RedstarSingleParticleMesonBlock; 
      int hsz = H.size(); 
      int psz = mom.size(); 
      inputList.resize(hsz * psz);   

      for (int h = 0; h < hsz; ++h)
        for(int p = 0; p < psz; ++p)
        {
          RedstarSingleParticleMesonInput *tmp; 
          tmp = new RedstarSingleParticleMesonInput;
          tmp->J = J; 
          tmp->H = H[h]; 
          tmp->parity = parity; 
          tmp->mom = mom[p]; 
          tmp->twoI_z = twoI_z; 
          tmp->name = name; 
          tmp->creation_op = creation_op; 
          tmp->smearedP = smearedP; 
          tmp->isProjected = isProjected; 
          tmp->t_slice = t_slice; 

          inputList[h*psz + p] = tmp;  
        }
    }

  //! xml writer
  void 
    RedstarSingleParticleMesonXML::write(ADATXML::XMLWriter &xml, const std::string &path ) const
    {
      ADATXML::push(xml,path);
      ADATXML::write(xml,"J",J);
      ADATXML::write(xml,"H",H);
      ADATXML::write(xml,"parity",parity);
      ADATXML::write(xml,"mom",mom);
      ADATXML::write(xml,"twoI_z",twoI_z);
      ADATXML::write(xml,"name",name);
      ADATXML::write(xml,"creation_op",creation_op);
      ADATXML::write(xml,"smearedP",smearedP);
      ADATXML::write(xml,"isProjected",isProjected); 
      ADATXML::write(xml,"t_slice",t_slice);
      ADATXML::pop(xml);
    }

} // radmat



