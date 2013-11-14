/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : redstar_unimproved_vector_current.cc

 * Purpose :

 * Creation Date : 11-11-2013

 * Last Modified : Thu 14 Nov 2013 03:14:27 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "redstar_unimproved_vector_current.h"
#include "redstar_single_particle_meson_block.h"
#include "radmat/construct_data/invert_subduction.h"

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
#include "radmat/utils/pow2assert.h"

#include "itpp/itbase.h"

#include <map>
#include <exception>
#include <sstream>

#define DEBUG_MSG_ON
#include "debug_props.h"


using namespace ADATXML; 

namespace radmat
{

  // a bunch of overloads and functions to move to cartesian coordinates

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
      doPrint(const ADATXML::Array<RedstarUnimprovedVectorCurrentPFrag> &t)
      {
        std::stringstream ss; 
        for(int i = 0; i < t.size(); ++i)
          ss << " " << t[i].coeff << " * " << t[i].name;
        return ss.str();  
      } 

    std::string 
      doPrint(const std::vector<RedstarUnimprovedVectorCurrentPFrag> &t)
      {
        std::stringstream ss; 
        for(int i = 0; i < t.size(); ++i)
          ss << " " << t[i].coeff << " * " << t[i].name;
        return ss.str();  
      } 

    // a debugging method
    void screen_dump (const EnsemRedstarBlock &k)
    {

#if 0
      std::cout << __func__ << ": nele = " << k.m_expr.size() << std::endl;

      //     for(int i = 0; i < k.m_expr.size(); ++i)
      //       std::cout << ensemFileName( k.m_expr[i].m_obj ) << std::endl;

      EnsemRedstarBlock::const_iterator it;
      for(it = k.begin(); it != k.end(); ++it)
      {
        std::cout << it->m_obj << std::endl;
        std::cout << ensemFileName(it->m_obj) << std::endl;
      }

#endif
    } 


    EnsemRedstarBlock operator*(const std::complex<double> &c, const EnsemRedstarBlock &l)
    {
      return SEMBLE::toScalar(c)*l; 
    }

    itpp::Vec<EnsemRedstarBlock> operator*(const itpp::Mat<std::complex<double> > &m,
        const itpp::Vec<EnsemRedstarBlock> &v)
    {
      POW2_ASSERT(v.size() == m.cols()); 
      itpp::Vec<EnsemRedstarBlock> ret(m.rows()); 

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

    EnsemRedstarBlock
      handle_temporal_work(const RedstarUnimprovedVectorCurrentInput &v)
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
        back.creation_op = v.creation_op; 
        back.smearedP = v.smearedP;
        back.isProjected = false; 
        back.t_slice = v.t_slice; 


        // loop over any photon flavors
        RedstarUnimprovedVectorCurrentInput::const_iterator it; 
        for(it = v.begin(); it != v.end(); ++it)
        {
          std::complex<double> weight(it->coeff,0.);
          back.name = it->name; 
          ret = ret + SEMBLE::toScalar(weight) * (piggy(&back)); 
        }

        DEBUG_MSG(exiting);

        return ret;
      }


    EnsemRedstarBlock
      handle_spatial_work(const RedstarUnimprovedVectorCurrentInput &v)
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
        back.creation_op = v.creation_op; 
        back.smearedP = v.smearedP;
        back.isProjected = false; 
        back.t_slice = v.t_slice; 

        // the transformation that takes us from helicity to cartesian coords
        itpp::Mat<std::complex<double> > M = invert2Cart(v.mom,v.creation_op); 
        itpp::Vec<EnsemRedstarBlock> cart, hel(3); 

        // loop over any photon flavors
        RedstarUnimprovedVectorCurrentInput::const_iterator it; 
        for(it = v.begin(); it != v.end(); ++it)
        {
          std::complex<double> weight(it->coeff,0.);
          back.name = it->name; 

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



    EnsemRedstarBlock 
      handle_work(const RedstarUnimprovedVectorCurrentInput &v)
      {
        if ( v.lorentz == 4 )
          return handle_temporal_work(v);   
        else if ( (v.lorentz > 0 ) && (v.lorentz < 4) )
          return handle_spatial_work(v); 
        else
        {
          std::cout << __PRETTY_FUNCTION__ 
            << ": Error, euclidean lorentz index is out of range (" 
            << v.lorentz << "), should be [1,4]" << std::endl; 
          exit(1); 
        }
      }



  } // anonomyous 


  std::string RedstarUnimprovedVectorCurrentInput::write(void) const
  {
    std::stringstream ss; 
    ss << "lorentz= " << lorentz << " mom= " << doPrint( mom ) 
      << " photons= " << doPrint( photons ) << " creation_op= " << creation_op 
      << " smearedP= " << smearedP << " t_slice= " << t_slice; 
    return ss.str(); 
  }


  EnsemRedstarBlock
    RedstarUnimprovedVectorCurrentBlock::operator()(const AbsRedstarInput_t *base) const
    {
      RedstarUnimprovedVectorCurrentInput dummy;
      POW2_ASSERT(base->type() == dummy.type());

      try
      {
        const RedstarUnimprovedVectorCurrentInput *derived;
        derived = dynamic_cast<const RedstarUnimprovedVectorCurrentInput*>(base); 
        dummy = *derived; 
      }
      catch(std::exception &e)
      {
        std::cout <<__PRETTY_FUNCTION__ << e.what() << std::endl; 
      }

      return handle_work(dummy); 
    }


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
        RedstarUnimprovedVectorCurrentXML::insertion &i,
        const int t_slice,
        const bool is_temporal)
    {
      RedstarUnimprovedVectorCurrentInput tmp; 
      int psz = i.photons.size(); 
      tmp.photons.resize(psz); 
      for(int idx = 0; idx < psz; ++idx)
        tmp.photons[idx] = i.photons[idx]; 
      tmp.creation_op = i.creation_op; 
      tmp.smearedP = i.smearedP; 
      tmp.t_slice = t_slice;

      if ( is_temporal ) 
      {
        RedstarUnimprovedVectorCurrentInput *t = new RedstarUnimprovedVectorCurrentInput;
        *t = tmp; 
        t->lorentz = 4; 
        v.push_back(rHandle<AbsRedstarInput_t>(t)); 
      }
      else
      {
        for (int lor = 1; lor < 4; ++lor)
        {
          RedstarUnimprovedVectorCurrentInput *t = new RedstarUnimprovedVectorCurrentInput;
          *t = tmp; 
          t->lorentz = lor; 
          v.push_back(rHandle<AbsRedstarInput_t>(t)); 
        }
      }

    }

    std::string toString(const RedstarUnimprovedVectorCurrentXML::insertion &i)
    {
      std::stringstream ss; 
      ss << "active= " << i.active << " create= " << i.creation_op 
        << " smear= " << i.smearedP << " photons: ";
      for(int j = 0; j < i.photons.size(); ++j)
        ss << i.photons[j].coeff << "x" << i.photons[j].name << "   ";
      return ss.str(); 
    }



  } // anonomyous  





  void write(ADATXML::XMLWriter &xml, 
      const std::string &path,
      const RedstarUnimprovedVectorCurrentPFrag &p)
  {
    ADATXML::push(xml,path);
    ADATXML::write(xml,"coeff",p.coeff);
    ADATXML::write(xml,"name",p.name); 
    ADATXML::pop(xml);
  }

  void write(ADATXML::XMLWriter &xml,
      const std::string &path, 
      const RedstarUnimprovedVectorCurrentXML::insertion &i)
  {
    ADATXML::push(xml,path);
    ADATXML::write(xml,"active",i.active);
    ADATXML::write(xml,"creation_op",i.creation_op);
    ADATXML::write(xml,"smearedP",i.smearedP); 
    write(xml,"photons",i.photons); 
    ADATXML::pop(xml);
  }

  void read(ADATXML::XMLReader &xml, 
      const std::string &path,
      RedstarUnimprovedVectorCurrentPFrag &p)  
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"coeff",p.coeff,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"name",p.name,__PRETTY_FUNCTION__); 
  }

  void read(ADATXML::XMLReader &xml, 
      const std::string &path,
      RedstarUnimprovedVectorCurrentXML::insertion &i)  
  {
    ADATXML::XMLReader ptop(xml,path); 
    doXMLRead(ptop,"active",i.active,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"creation_op",i.creation_op,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"smearedP",i.smearedP,__PRETTY_FUNCTION__); 
    doXMLRead(ptop,"photons",i.photons,__PRETTY_FUNCTION__); 
  }

  void 
    RedstarUnimprovedVectorCurrentXML::read(ADATXML::XMLReader &xml, const std::string &path)
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"pmin",pmin,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"pmax",pmax,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"t_slice",t_slice,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"time",time,__PRETTY_FUNCTION__);       
      doXMLRead(ptop,"space",space,__PRETTY_FUNCTION__);       


      objFunctorPtr = rHandle<AbsRedstarBlock_t>(new RedstarUnimprovedVectorCurrentBlock); 

      // populate the inputList vector
      read_in_photons(inputList,time,t_slice,true); 
      read_in_photons(inputList,space,t_slice,false);  
    }

  std::string 
    RedstarUnimprovedVectorCurrentXML::write(void) const
    {
      std::stringstream ss;
      ss << "pmin= " << pmin << " pmax=" << pmax << " t_slice= " << t_slice; 
      ss << "\ntime:\n" << toString(time) << std::endl;
      ss << "\nspace:\n" << toString(space) << std::endl;
      return ss.str(); 
    }

  void 
    RedstarUnimprovedVectorCurrentXML::write(
        ADATXML::XMLWriter &xml, 
        const std::string &path) const
    {
      ADATXML::push(xml,path);
      ADATXML::write(xml,"pmin",pmin);
      ADATXML::write(xml,"pmax",pmax);
      ADATXML::write(xml,"t_slice",t_slice);
      ::radmat::write(xml,"time",time);
      ::radmat::write(xml,"space",space); 
      ADATXML::pop(xml); 
    }

} // radmat
