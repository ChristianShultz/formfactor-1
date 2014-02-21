#ifndef LORENTZFF_PIPI_H_H_GUARD
#define LORENTZFF_PIPI_H_H_GUARD


#include "lorentzff_polarization_embedding.h"
#include "radmat/utils/pow2assert.h"
#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_canonical_frame_formfactor.h"
#include "radmat/utils/stringify.h"
#include <complex>

namespace radmat
{


  // only one ff
  struct PiPiF1impl
    : public FormFacRotationManager<PiPiF1impl, std::complex<double> >
  {
    virtual ~PiPiF1impl() {}

    virtual  std::string ff_impl(void) const
    {
      return std::string(" F_1(Q^2) p_+^{\\mu} ");
    }

    // return a complex version of p_+
    virtual Tensor<std::complex<double> , 1> 
      impl(const Tensor<double,1> &p_f,
          const Tensor<double,1> &p_i,
          const double mom_fac,
          int h_f,
          int h_i) const
      {
        return convertTensorUnderlyingType<std::complex<double>,double,1>( pPlus(p_f,p_i) );
      }
  };

  struct PiPiF1; 
  REGISTER_STRINGIFY_TYPE( PiPiF1 ); 

  struct PiPiF1
    : public canonicalFrameFormFactor<0,0,PiPiF1impl>
  {
    virtual ~PiPiF1() {} 
    virtual std::string id() const {return Stringify<PiPiF1>();}
  };



  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////


  // generate a list for the PiPi constructor
  template<int embedl, int embedr>
    FFAbsBase_t::FFAbs_list PiPiGenList(void)
    {
      FFAbsBase_t::FFAbs_list retPiPi;
      FFAbsBase_t::BBType *blockPtr;
      blockPtr = new PiPiF1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPi.push_back(FFAbsBase_t::BBHandle_t(blockPtr));
      return retPiPi;
    }


  ////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////

  // can only instantiate the <0,0> type now
  template<int embedl, int embedr> struct PiPi; 
  REGISTER_STRINGIFY_TYPE2( PiPi<0,0> ); 


  // only need to derive the constructor, everything else is in the 
  // base class, do this polymorphically, make some function that 
  // returns the appropriate handle based on the requested matrix element type
  template<int embedl, int embedr>
    struct PiPi : public FFAbsBase_t
  {
    PiPi(void)
      : FFAbsBase_t(radmat::PiPiGenList<embedl,embedr>())  
    {  } 

    PiPi& operator=(const PiPi &o)
    {

      if(this != &o)
        FFAbsBase_t::operator=(o);
      return *this;
    }

    // no slicing
    PiPi(const PiPi &o)
      : FFAbsBase_t(o)
    {  }

    virtual ~PiPi() {} 

    virtual std::string id(void) const { return Stringify< PiPi<embedl,embedr> >(); }
    virtual int left_spin(void) const {return embedl;}
    virtual int right_spin(void) const {return embedr;}


    private:
    // I'm not sure if these could inherit so we will hide them as well
    PiPi(const FFAbsBase_t::FFAbs_list &);
    PiPi(const FFAbsBase_t::FFAbs_list);
  };


} // close radmat

#endif
