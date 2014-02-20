#ifndef LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H
#define LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H 


#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"

namespace radmat
{

  // this confusing thing basically hides the spin from the 
  // form factor classes -- too much effort to recode it for
  // easier reading 
  template<int Jl, int Jr, class DerivedFF> 
    struct canonicalFrameFormFactor
    : public DerivedFF,
    public FFAbsBlockBase_t
  {
    typedef MomHelPair_t p4_t; 

    virtual ~canonicalFrameFormFactor() {}

    virtual Tensor<std::complex<double>,1>
      operator()(const p4_t &l, const p4_t &r, const double kick) const 
      {
        return DerivedFF::operator()(l.first,r.first,kick,Jl,Jr,l.second,r.second); 
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
  //>
  // struct F1 
  // : public canonicalFrameFormFactor<1,1,F1impl> 
  // {
  //  virtual ~F1() {}
  // };
  //
  /////////////////////////////////////////////////////////////////////////

} // radmat 

#endif /* LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H */
