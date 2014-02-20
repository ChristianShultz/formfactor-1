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

  namespace PiPi
  {

    // only one ff
    struct F1impl
      : public FormFacRotationManager<F1impl,<std::complex<double> >
    {
      virtual ~F1impl() {}


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


    struct F1
      : public canonicalFrameFormFactor<0,0,F1impl>
    {
      virtual ~F1() {} 
    };



    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////


    // generate a list for the PiPi constructor
    template<int embed>
    FFAbsBase_t::FFAbs_list genList(void)
    {
      FFAbsBase_t::FFAbs_list retPiPi;
      FFAbsBase_t::BBType *blockPtr;
      blockPtr = new F1();
      POW2_ASSERT(blockPtr);  // blow up if something went wrong
      retPiPi.push_back(FFAbsBase_t::BBHandle_t(blockPtr));
      return retPiPi;
    }



    ////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////

    // embed spin into operator name -- only 0,0 compiles
    template<int embedl, int embedr> struct PiPi; 
    REGISTER_STRINGIFY_TYPE( PiPi<0,0> ); 


    // only need to derive the constructor, everything else is in the 
    // base class, do this polymorphically, make some function that 
    // returns the appropriate handle based on the requested matrix element type
    template<int embedl, int embedr>
    struct PiPi : public FFAbsBase_t
    {
      PiPi(void)
        : FFAbsBase_t(radmat::PiPi::genList<embed>())  
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
      virtual int left_spin const {return embedl;}
      virtual int right_spin const {return embedr;}


      private:
      // I'm not sure if these could inherit so we will hide them as well
      PiPi(const FFAbsBase_t::FFAbs_list &);
      PiPi(const FFAbsBase_t::FFAbs_list);
    };


  } // close PiPi

} // close radmat

#endif
