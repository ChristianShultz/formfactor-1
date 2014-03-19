/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_improved_vector_current.cc

 * Purpose :

 * Creation Date : 20-11-2013

 * Last Modified : Tue 18 Mar 2014 09:25:40 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "redstar_improved_vector_current.h"
#include "redstar_cartesian_interface.h"
#include "redstar_single_particle_meson_block.h"
#include "redstar_invert_subduction.h"
#include "redstar_photon_props.h"
#include "semble/semble_meta.h"
#include "adat/map_obj.h"
#include "ensem/ensem.h"
#include "hadron/ensem_filenames.h"
#include "hadron/irrep_util.h"
#include "hadron/clebsch.h"
#include "formfac/formfac_qsq.h"
#include "redstar_photon_polarization_tensor.h"
#include "radmat/utils/tensor.h"
#include "radmat/utils/obj_expr_t.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/mink_qsq.h"
#include "itpp/itbase.h"
#include <map>
#include <exception>
#include <sstream>

#define DEBUG_MSG_OFF
#include "radmat/utils/debug_handler.h"


using namespace ADATXML; 


namespace radmat
{

  namespace
  {

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

    std::string 
      doPrint(const ADATXML::Array<RedstarImprovedVectorCurrentPFrag> &t)
      {
        std::stringstream ss; 
        int sz = t.size(); 
        for(int i = 0; i < sz; ++i)
          ss << " " << t[i].op_coeff << " * " << t[i].op_name
            << "\n       + <"<< t[i].i_coeff_r << ","
            << t[i].i_coeff_i << "> * " << t[i].i_name 
            << "\nm_qxa_t= "<<  t[i].m_qxa_t;
        return ss.str();  
      } 

    std::string 
      doPrint(const std::vector<RedstarImprovedVectorCurrentPFrag> &t)
      {
        std::stringstream ss; 
        int sz = t.size(); 
        for(int i = 0; i < sz; ++i)
          ss << " " << t[i].op_coeff << " * " << t[i].op_name
            << "\n       + <"<< t[i].i_coeff_r << ","
            << t[i].i_coeff_i << "> * " << t[i].i_name ;
        return ss.str();  
      } 

    std::string 
      doPrint(const RedstarImprovedVectorCurrentParPack &p)
      {
        std::stringstream ss; 
        ss << "xi= " << p.xi << " L_s= " << p.L_s
          << " d_1= " << p.d_1 << " gamma_f= " << p.gamma_f 
          << " mu_tilde_t= " << p.mu_tilde_t
          << " msnk= " << p.msnk << " msrc= " << p.msrc;
        return ss.str(); 
      }

    std::string 
      doPrint(const RedstarImprovedVectorCurrentInsertion &f)
      {
        std::stringstream ss;
        ss << "active= " << f.active << " creation_op= " << PHOTON_CREATE
          << " smearedP= " << f.smearedP << "\nphotons:\n" << doPrint(f.photons); 
        return ss.str(); 
      }

    // a debugging method
    void screen_dump (const EnsemRedstarBlock &k)
    {
#if 1 
      std::cout << __func__ << ": nele = " << k.m_expr.size() << std::endl;

      //     for(int i = 0; i < k.m_expr.size(); ++i)
      //       std::cout << ensemFileName( k.m_expr[i].m_obj ) << std::endl;

      EnsemRedstarBlock::const_iterator it;
      for(it = k.begin(); it != k.end(); ++it)
      {
        std::cout << it->m_obj << std::endl;
        // std::cout << Hadron::ensemFileName(it->m_obj) << std::endl;
      }

#endif
    } 


    std::string stringy_mom(const ADATXML::Array<int> mom)
    {
      ADATXML::Array<int> can = FF::canonicalOrder(mom);
      std::stringstream ss;
      ss << "p" << can[0] << can[1] << can[2];
      return ss.str(); 
    }


  } // anonomyous 



  //////////////////////////////////////////////////////
  //                INPUT STUFF                       //
  //////////////////////////////////////////////////////

  //////////////////////////////////////////////////////
  //
  std::string
    RedstarImprovedVectorCurrentInput::write(void) const
    {
      std::stringstream ss; 
      ss << "lorentz= " << lorentz << " mom= " << doPrint( mom ) 
        << " psrc= " << doPrint(psrc) << " psnk= " << doPrint(psnk)
        << " photons= " << doPrint( photons ) << " creation_op= " << PHOTON_CREATE 
        << " smearedP= " << smearedP << " t_slice= " << t_slice
        << "\nipack:\n" << doPrint(ipack);
      return ss.str(); 

    }


  //////////////////////////////////////////////////////
  //                BLOCK STUFF                       //
  //////////////////////////////////////////////////////

  namespace
  {

    //////////////////////////////////////////////////////
    //
    itpp::Vec<EnsemRedstarBlock>
      cartesianizeImprovement(const RedstarImprovedVectorCurrentInput &v)
      {
        DEBUG_MSG(entering);
        itpp::Vec<EnsemRedstarBlock> ret(3); 

        // piggy back off single particle mesons
        RedstarSingleParticleMesonInput back; 
        RedstarSingleParticleMesonBlock piggy;
        back.J = 1; 
        back.parity = false; 
        back.mom = v.mom; 
        back.twoI_z = 0; 
        back.creation_op = PHOTON_CREATE; 
        back.smearedP = v.smearedP;
        back.isProjected = false; 
        back.t_slice = v.t_slice; 

        // the transformation that takes us from helicity to cartesian coords
        itpp::Mat<std::complex<double> > M = invert2Cart(v.mom,back.creation_op); 
        itpp::Vec<EnsemRedstarBlock> cart, hel(3); 

        // loop over any photon flavors
        RedstarImprovedVectorCurrentInput::const_iterator it; 
        for(it = v.begin(); it != v.end(); ++it)
        {

          if(it->i_coeff_r == 0.)
            if(it->i_coeff_i == 0.)
              continue;  

          // paranoia on my part
          if( it->op_coeff == 0.)
          {
            std::cout << __PRETTY_FUNCTION__ 
              << ": stupidity -- what are you doing improving a photon with"
              << " a coefficient of zero -- are you sure this" 
              << " is what you intend? " << std::endl; 
          }

          std::complex<double> weight(it->i_coeff_r,it->i_coeff_i);
          back.name = it->i_name; 

          // little loop over helicity 
          for (int h = -1; h < 2; ++h)  
          {
            back.H = h; 
            hel[1 -h] = piggy(&back); 
          }

          // transform to cartesian coordinates
          ret = ret + weight * ( M * hel );
        }

        DEBUG_MSG(exiting); 

       //  for(int i = 0; i < ret.size(); ++i)
       //   screen_dump(ret[i]); 

        return ret; 
      }




    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock
      handle_improved_temporal_work(const RedstarImprovedVectorCurrentInput &v)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        POW2_ASSERT( v.lorentz == 4 );

        itpp::Vec<EnsemRedstarBlock> cart = cartesianizeImprovement(v); 

        for(int idx = 0; idx < 3; ++idx)
        {
          if ( v.mom[idx] == 0 ) 
            continue; 

          if ( cart[idx].begin() == cart[idx].end() ) 
            continue; 
          // a_s * d_1 * 2pi * p_j  
          //    = d_1 * 2pi * n_j / L_s
          double weight =  v.ipack.d_1 * 2. * acos(-1) * double(v.mom[idx]) / double(v.ipack.L_s);

          // only overloaded mul for std::complex<double> -- 
          //    just hack it
          ret = ret + std::complex<double>(weight,0.) * cart[idx];
        }

        DEBUG_MSG(exiting); 
        return ret; 
      }



    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock
      handle_standard_temporal_work(const RedstarImprovedVectorCurrentInput &v)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        POW2_ASSERT( v.lorentz == 4 );

        // piggy back off single particle mesons
        RedstarSingleParticleMesonInput back; 
        RedstarSingleParticleMesonBlock piggy;
        back.J = 0; 
        back.H = 0; 
        back.parity = true; 
        back.mom = v.mom; 
        back.twoI_z = 0; 
        back.creation_op = PHOTON_CREATE; 
        back.smearedP = v.smearedP;
        back.isProjected = false; 
        back.t_slice = v.t_slice; 


        // loop over any photon flavors
        RedstarImprovedVectorCurrentInput::const_iterator it; 
        for(it = v.begin(); it != v.end(); ++it)
        {
          std::complex<double> weight(it->op_coeff,0.);
          back.name = it->op_name; 
          ret = ret + SEMBLE::toScalar(weight) * (piggy(&back)); 
        }

        DEBUG_MSG(exiting);

        return ret;
      }
#if 0
    // This choice corresponds to 
    //      g_k -> g_k * ( 1 - 2* d_1*gamma_f * mu_tilde * atmq ) 
    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock
      handle_standard_spatial_work(const RedstarImprovedVectorCurrentInput &v)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        POW2_ASSERT( v.lorentz < 4 );
        POW2_ASSERT( v.lorentz > 0 );

        double standard_i_coeff; 
        standard_i_coeff = 2.*v.ipack.d_1*v.ipack.gamma_f*v.ipack.mu_tilde_t; 

        // piggy back off single particle mesons
        RedstarSingleParticleMesonInput back; 
        RedstarSingleParticleMesonBlock piggy;
        back.J = 1; 
        back.parity = false; 
        back.mom = v.mom; 
        back.twoI_z = 0; 
        back.creation_op = v.creation_op; 
        back.smearedP = v.smearedP;
        back.isProjected = false; 
        back.t_slice = v.t_slice; 

        // the transformation that takes us from helicity to cartesian coords
        itpp::Mat<std::complex<double> > M = invert2Cart(v.mom,v.creation_op); 
        itpp::Vec<EnsemRedstarBlock> cart, hel(3); 

        // loop over any photon flavors
        RedstarImprovedVectorCurrentInput::const_iterator it; 
        for(it = v.begin(); it != v.end(); ++it)
        {
          // the isospin weighting (based on quark charges and 
          //     operator definitions)
          std::complex<double> weight(it->op_coeff,0.);
          // the weight from the improvement term 
          //  gamma_mu -> gamma_mu * ( 1 - 2*a_s*d_1*mu*m_0*gamma_f/xi)
          weight *= ( 1. - it->m_qxa_t * standard_i_coeff);  
          back.name = it->op_name; 

          // little loop over helicity 
          for (int h = -1; h < 2; ++h)  
          {
            back.H = h; 
            hel[1 -h] = piggy(&back); 
          }

          // transform to cartesian coordinates
          cart = M * hel;

          // NB: some numbering assumptions in here -- pull the guy we want
          ret = ret + SEMBLE::toScalar(weight) * ( cart[ v.lorentz - 1] ) ;   
        }

        DEBUG_MSG(exiting);

        return ret;
      }

    double 
      performBoost(const double m, 
          const ADATXML::Array<int> mom, 
          const double xi, 
          const int L_s)
      {
        double unit,p;
        for(int i = 0; i < 3; ++i)
          p += double(mom[i]*mom[i]); 

        unit = mom_factor(xi,L_s); 
        p*= ( unit * unit );

        return sqrt( m*m + p ); 
      }
#endif 

    // This choice corresponds to absorbing that extra factor
    // into the non-perturbative Z_V
    //
    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock
      handle_standard_spatial_work(const RedstarImprovedVectorCurrentInput &v)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        POW2_ASSERT( v.lorentz < 4 );
        POW2_ASSERT( v.lorentz > 0 );

        // piggy back off single particle mesons
        RedstarSingleParticleMesonInput back; 
        RedstarSingleParticleMesonBlock piggy;
        back.J = 1; 
        back.parity = false; 
        back.mom = v.mom; 
        back.twoI_z = 0; 
        back.creation_op = PHOTON_CREATE; 
        back.smearedP = v.smearedP;
        back.isProjected = false; 
        back.t_slice = v.t_slice; 

        // the transformation that takes us from helicity to cartesian coords
        itpp::Mat<std::complex<double> > M = invert2Cart(v.mom,back.creation_op); 
        itpp::Vec<EnsemRedstarBlock> cart, hel(3); 

        // loop over any photon flavors
        RedstarImprovedVectorCurrentInput::const_iterator it; 
        for(it = v.begin(); it != v.end(); ++it)
        {
          // the isospin weighting (based on quark charges and 
          //     operator definitions)
          std::complex<double> weight(it->op_coeff,0.);
          back.name = it->op_name; 

          // little loop over helicity 
          for (int h = -1; h < 2; ++h)  
          {
            back.H = h; 
            hel[1 -h] = piggy(&back); 
          }

          // transform to cartesian coordinates
          cart = M * hel;

          // NB: some numbering assumptions in here -- pull the guy we want
          ret = ret + SEMBLE::toScalar(weight) * ( cart[ v.lorentz - 1] ) ;   
        }

        DEBUG_MSG(exiting);

        return ret;
      }

    //////////////////////////////////////////////////////
    //
    double 
      performBoost(const double m, 
          const ADATXML::Array<int> mom, 
          const double xi, 
          const int L_s)
      {
        double unit,p;
        for(int i = 0; i < 3; ++i)
          p += double(mom[i]*mom[i]); 

        unit = mom_factor(xi,L_s); 
        p*= ( unit * unit );

        return sqrt( m*m + p ); 
      }


    // Follows the conventions of the momentumTransfer 
    //          function, pretty sure this gets the 
    //          sign correct but need to think 
    //          about it some more before I believe it
    //
    //////////////////////////////////////////////////////
    //
    double  
      deltaE(const RedstarImprovedVectorCurrentInput &i)
      { 
        double Esnk = performBoost(i.ipack.msnk,i.psnk,i.ipack.xi,i.ipack.L_s);  
        double Esrc = performBoost(i.ipack.msrc,i.psrc,i.ipack.xi,i.ipack.L_s);  

        double DE = Esrc - Esnk; 

        if ( i.creation_op ) 
          return -DE; 

        return DE; 
      }


    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock
      handle_improved_spatial_work(const RedstarImprovedVectorCurrentInput &v)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        POW2_ASSERT( v.lorentz < 4 );
        POW2_ASSERT( v.lorentz > 0 );

        itpp::Vec<EnsemRedstarBlock> cart = cartesianizeImprovement(v); 

        double dw = v.ipack.d_1 * v.ipack.gamma_f; 
        dw *= deltaE(v);

        ret = std::complex<double>(dw ,0.) * cart[v.lorentz -1];

        DEBUG_MSG(exiting); 
        return ret; 
      }

    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock 
      handle_temporal_work(const RedstarImprovedVectorCurrentInput &i)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        ret =  handle_standard_temporal_work(i)
          + handle_improved_temporal_work(i); 

        DEBUG_MSG(exiting);

        return ret; 
      }

    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock 
      handle_spatial_work(const RedstarImprovedVectorCurrentInput &i)
      {
        DEBUG_MSG(entering);
        EnsemRedstarBlock ret; 

        ret = handle_standard_spatial_work(i)
          + handle_improved_spatial_work(i); 

        DEBUG_MSG(exiting);

        return ret; 
      }

    //////////////////////////////////////////////////////
    //
    EnsemRedstarBlock 
      handle_work(const RedstarImprovedVectorCurrentInput &i)
      {
        if ( i.lorentz == 4 ) 
          return handle_temporal_work(i); 
        else 
          return handle_spatial_work(i); 
      }
  } // anonomyous 

  //////////////////////////////////////////////////////
  //
  EnsemRedstarBlock
    RedstarImprovedVectorCurrentBlock::operator()(const AbsRedstarInput_t *base) const
    {
      RedstarImprovedVectorCurrentInput dummy; 
      POW2_ASSERT(base->type() == dummy.type()); 

      try
      {
        const RedstarImprovedVectorCurrentInput *derived; 
        derived = dynamic_cast<const RedstarImprovedVectorCurrentInput*>(base); 
        dummy = *derived; 
      }
      catch (std::exception &e)
      {
        std::cout << __PRETTY_FUNCTION__ << e.what() << std::endl;
      }

      return handle_work(dummy); 
    } 



  //////////////////////////////////////////////////////
  //                 XML STUFF                        //
  //////////////////////////////////////////////////////

  namespace
  {
    template<typename T>
      void doXMLRead(ADATXML::XMLReader &ptop, 
          const std::string &path, 
          T &place, 
          const char * f)
      {
        if(ptop.count(path) > 0)
          read(ptop,path,place);
        else
        {
          std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " 
            << f << " trying to read path, " << path
            << ", path was empty, exiting" << std::endl;
          exit(1);
        }
      }


    void read_in_photons(std::vector<rHandle<AbsRedstarInput_t> > &v,
        const RedstarImprovedVectorCurrentInsertion &i, 
        const RedstarImprovedVectorCurrentParPack & ipack,
        const int t_slice, 
        const bool is_temporal)
    {
      RedstarImprovedVectorCurrentInput tmp; 
      int psz = i.photons.size(); 
      tmp.photons.resize(psz); 
      for(int idx = 0; idx < psz; ++idx)
        tmp.photons[idx] = i.photons[idx];
      tmp.creation_op = i.creation_op; 
      tmp.smearedP = i.smearedP; 
      tmp.t_slice = t_slice; 
      tmp.ipack = ipack; 


      if ( is_temporal ) 
      {
        RedstarImprovedVectorCurrentInput *t;
        t = new RedstarImprovedVectorCurrentInput(tmp); 
        t->lorentz = 4;
        v.push_back(rHandle<AbsRedstarInput_t>(t)); 
      }
      else
      {
        for (int lor = 1 ; lor < 4; ++lor)
        {
          RedstarImprovedVectorCurrentInput *t;
          t = new RedstarImprovedVectorCurrentInput(tmp); 
          t->lorentz = lor;
          v.push_back(rHandle<AbsRedstarInput_t>(t)); 
        }
      }

    }






  } // anonomyous 

  //////////////////////////////////////////////////////
  //
  void 
    RedstarImprovedVectorCurrentXML::read(ADATXML::XMLReader &xml, 
        const std::string &path)
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"pmin",pmin,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"pmax",pmax,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"t_slice",t_slice,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"ipack",ipack,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"time",time,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"space",space,__PRETTY_FUNCTION__); 


      objFunctorPtr = rHandle<AbsRedstarBlock_t>(new RedstarImprovedVectorCurrentBlock); 

      if ( time.active ) 
        read_in_photons(inputList,time,ipack,t_slice,true); 
      if ( space.active )
        read_in_photons(inputList,space,ipack,t_slice,false);

      if( time.active && space.active )
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ 
          << ": Warning, you decided to mix time and"
          << " space and I don't think you should" << std::endl; 
      }


      if ( inputList.empty() ) 
      {
        std::cout << __PRETTY_FUNCTION__ << __FILE__ << __LINE__
          << ": Warning, no photons were read" << std::endl;
      }
    }

  std::string 
    RedstarImprovedVectorCurrentXML::write(void) const
    {
      std::stringstream ss; 
      ss << "pmin= " << pmin << " pmax= " << pmax 
        << " t_slice= " << t_slice << "\nipack: " << doPrint(ipack)  
        << "\ntime: " << doPrint(time) << "\nspace: " << doPrint(space); 
      return ss.str(); 
    }

  void 
    RedstarImprovedVectorCurrentXML::write(ADATXML::XMLWriter &xml, 
        const std::string &path) const
    {
      ADATXML::push(xml,path); 
      ADATXML::write(xml,"pmin",pmin);
      ADATXML::write(xml,"pmax",pmax);
      ADATXML::write(xml,"t_slice",t_slice); 
      ::radmat::write(xml,"ipack",ipack);
      ::radmat::write(xml,"time",time);
      ::radmat::write(xml,"space",space);
      ADATXML::pop(xml); 
    } 

  ///////////////////////////////////////////////////////////////////////////
  //
  //  XML 

  ///////////////////////////////////////////////////////////////////////////
  // xml readers
  void read(ADATXML::XMLReader &xml,
      const std::string &path,
      RedstarImprovedVectorCurrentPFrag &f)
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"op_name",f.op_name,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"op_coeff",f.op_coeff,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"i_name",f.i_name,__PRETTY_FUNCTION__);

    // break condition
    if ( f.i_name == "no_improvement" )
    {
      f.i_coeff_r = 1000000.;
      f.i_coeff_i = 1000000.;
      f.m_qxa_t = 0.;
      return;
    }
    doXMLRead(ptop,"i_coeff_r",f.i_coeff_r,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"i_coeff_i",f.i_coeff_i,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"m_qxa_t",f.m_qxa_t,__PRETTY_FUNCTION__); 
  }

  ///////////////////////////////////////////////////////////////////////////
  //
  void read(ADATXML::XMLReader &xml,
      const std::string &path,
      RedstarImprovedVectorCurrentInsertion &f)
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"active",f.active,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"smearedP",f.smearedP,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"photons",f.photons,__PRETTY_FUNCTION__);
  }

  ///////////////////////////////////////////////////////////////////////////
  //
  void read(ADATXML::XMLReader &xml,
      const std::string &path,
      RedstarImprovedVectorCurrentParPack &f)
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"xi",f.xi,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"L_s",f.L_s,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"d_1",f.d_1,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"gamma_f",f.gamma_f,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"mu_tilde_t",f.mu_tilde_t,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"msnk",f.msnk,__PRETTY_FUNCTION__);
    doXMLRead(ptop,"msrc",f.msrc,__PRETTY_FUNCTION__);
  }

  ///////////////////////////////////////////////////////////////////////////
  // xml writers
  void write(ADATXML::XMLWriter &xml,
      const std::string &path,
      const RedstarImprovedVectorCurrentPFrag &f)
  {
    ADATXML::push(xml,path);  
    ADATXML::write(xml,"op_name",f.op_name);
    ADATXML::write(xml,"op_coeff",f.op_coeff);
    ADATXML::write(xml,"i_name",f.i_name);
    ADATXML::write(xml,"i_coeff_r",f.i_coeff_r);
    ADATXML::write(xml,"i_coeff_i",f.i_coeff_i);
    ADATXML::write(xml,"m_qxa_t",f.m_qxa_t); 
    ADATXML::pop(xml);
  }

  ///////////////////////////////////////////////////////////////////////////
  //
  void write(ADATXML::XMLWriter &xml,
      const std::string &path,
      const RedstarImprovedVectorCurrentInsertion &f)
  {
    ADATXML::push(xml,path);
    ADATXML::write(xml,"active",f.active);
    ADATXML::write(xml,"smearedP",f.smearedP);
    ADATXML::write(xml,"photons",f.photons);
    ADATXML::pop(xml);
  }

  ///////////////////////////////////////////////////////////////////////////
  //
  void write(ADATXML::XMLWriter &xml,
      const std::string &path,
      const RedstarImprovedVectorCurrentParPack &f)
  {
    ADATXML::push(xml,path);
    ADATXML::write(xml,"xi",f.xi);
    ADATXML::write(xml,"L_s",f.L_s);
    ADATXML::write(xml,"d_1",f.d_1);
    ADATXML::write(xml,"gamma_f",f.gamma_f);
    ADATXML::write(xml,"mu_tilde_t",f.mu_tilde_t);
    ADATXML::write(xml,"msnk",f.msnk);
    ADATXML::write(xml,"msrc",f.msrc);
    ADATXML::pop(xml); 
  }

} // radmat
