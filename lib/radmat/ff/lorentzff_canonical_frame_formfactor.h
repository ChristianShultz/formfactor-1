#ifndef LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H
#define LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H 


#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"

namespace radmat
{

  template<int Jl, int Jr, int hl, int hr, class DerivedFF>
    struct canonicalFrameFormFactor
    : public DerivedFF,
    public ffBlockBase_t<std::complex<double> >
  {
    typedef Tensor<double,1> p4_t; 

    virtual ~embedRotationManager() {}

    virtual Tensor<std::complex<double>,1>
      operator()(const p4_t &l, const p4_t &r, const double kick) const 
      {
        return DerivedFF::operator()(l,r,kick,Jl,Jr,hl,hr); 
      }   

    virtual std::string 
      ff(void) const
      {
        return DerivedFF::ff_impl(); 
      }
  };


  /////////////////////////////////////////////////////////////////////////
  //
  //  struct F1impl
  //  : public FormFacRotationManager<F1impl>
  //  {
  //    virtual ~F1impl() {}
  //    virtual std::string ff_impl() const { do stuff}
  //    virtual Tensor impl(args) const { do other stuff}
  //  };
  //
  // template<int lambda_l, int lambda_r>
  // struct F1 
  // : public canonicalFrameFormFactor<1,1,lambda_l, lambda_r, RhoRhoImpl> 
  // {
  //  virtual ~F1() {}
  // };
  //
  /////////////////////////////////////////////////////////////////////////

} // radmat 

#endif /* LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H */
