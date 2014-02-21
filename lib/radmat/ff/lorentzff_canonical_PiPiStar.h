#ifndef LORENTZFF_PIPISTAR_H_H_GUARD
#define LORENTZFF_PIPISTAR_H_H_GUARD


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/pow2assert.h"
#include "radmat/utils/stringify.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_frame_formfactor.h"
#include <complex>

namespace radmat
{

  // only one ff
  struct PiPiStarF1impl
    : public FormFacRotationManager<PiPiStarF1impl, std::complex<double> >
  {
    virtual ~PiPiStarF1impl() {}

    virtual  std::string ff_impl(void) const
    {
      return std::string("F_1(Q^2)\\left(- p_+^{\\mu}\frac{Q^2}{m_{\\pi*}^2 -m_{\\pi}^2} + p_-\\right)");
    }

    // return a complex version of p_+
    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f,
          const Tensor<double,1> &p_i,
          const double mom_fac,
          int h_f,
          int h_i) const
      {
        Tensor<std::complex<double>,1> pp,pm;
        pp = convertTensorUnderlyingType<std::complex<double>,double,1>(pPlus(p_f,p_i));
        pm = convertTensorUnderlyingType<std::complex<double>,double,1>(pMinus(p_f,p_i));

        // double num = (p_f - p_i) * (g_dd() * (p_f - p_i));
        // double denom = (p_f - p_i) * (g_dd() * (p_f + p_i));

        double num = value(contract(p_f-p_i,p_f-p_i,g_dd(),0,0));
        double denom = value(contract(p_f-p_i,p_f+p_i,g_dd(),0,0));


        POW2_ASSERT_DEBUG(fabs(denom) > 1e-14);

        return ( -(num/denom)*pp +  pm);
      }
  };


  struct PiPiStarF1;
  REGISTER_STRINGIFY_TYPE( PiPiStarF1 ); 

  struct PiPiStarF1
    : public canonicalFrameFormFactor<0,0,PiPiStarF1impl>
  {
    virtual ~PiPiStarF1() {} 
    std::string id() const { return Stringify<PiPiStarF1>(); }
  };



  // generate a list for the PiPiStar constructor
  template<int embedl, int embedr>
    FFAbsBase_t::FFAbs_list PiPiStarGenList(void)
    {
      FFAbsBase_t::FFAbs_list retPiPiStar;
      FFAbsBase_t::BBType *blockPtr;
      blockPtr = new PiPiStarF1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPiStar.push_back(FFAbsBase_t::BBHandle_t(blockPtr));
      return retPiPiStar;
    }


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  template<int embedl, int embedr> struct PiPiStar; 
  REGISTER_STRINGIFY_TYPE2( PiPiStar<0,0> ); 

  // only need to derive the constructor, everything else is in the 
  // base class, do this polymorphically, make some function that 
  // returns the appropriate handle based on the requested matrix element type
  template<int embedl, int embedr>
    struct PiPiStar : public FFAbsBase_t
  {
    PiPiStar(void)
      : FFAbsBase_t(radmat::PiPiStarGenList<embedl,embedr>())  
    {  } 

    PiPiStar& operator=(const PiPiStar &o)
    {

      if(this != &o)
        FFAbsBase_t::operator=(o);
      return *this;
    }

    // no slicing
    PiPiStar(const PiPiStar &o)
      : FFAbsBase_t(o)
    {  }

    virtual std::string id(void) const { return Stringify< PiPiStar<embedl,embedr> >(); }
    virtual int left_spin(void) const {return embedl;}
    virtual int right_spin(void) const {return embedr;}

    private:
    // I'm not sure if these could inherit so we will hide them as well
    PiPiStar(const FFAbsBase_t::FFAbs_list &);
    PiPiStar(const FFAbsBase_t::FFAbs_list);
  };

} // close radmat

#endif /* PIPISTAR GUARD */
