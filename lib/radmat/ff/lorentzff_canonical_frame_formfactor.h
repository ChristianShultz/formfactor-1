#ifndef LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H
#define LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H 


#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"

namespace radmat
{

  template<int Jl, int Jr, int hl, int hr, class DerivedFF, bool is_ff>
    struct canonicalFrameFormFactor
    : public DerivedFF,
    public ffBlockBase_t<std::complex<double> >
  {
    typedef Tensor<double,1> p4_t; 

    virtual ~canonicalFrameFormFactor() {}

    virtual Tensor<std::complex<double>,1>
      operator()(const p4_t &l, const p4_t &r, const double kick) const 
      {
        //
        // fix up any lattice mass splittings in the case of a form factor 
        //    (not a transition )
        //
        if ( is_ff ) 
        {   
          Tensor<double,0> ml,mr; 
          double m; 
          ml = contract(l,contract(g_dd(),l,1,0),0,0);
          mr = contract(r,contract(g_dd(),r,1,0),0,0);
          m = 0.5*(ml.value() + mr.value()); 
          Tensor<double,1> ll,rr; 
          ll = l; 
          rr = r; 
          ll[0] = sqrt( m + l[1]*l[1] + l[2]*l[2] + l[3]*l[3]);
          rr[0] = sqrt( m + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);

          return DerivedFF::operator()(ll,rr,kick,Jl,Jr,hl,hr); 
        }

        // dont do anything for transitions 
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
  // : public canonicalFrameFormFactor<1,1,lambda_l, lambda_r, RhoRhoImpl, bool> 
  // {
  //  virtual ~F1() {}
  // };
  //
  /////////////////////////////////////////////////////////////////////////

} // radmat 

#endif /* LORENTZFF_CANONICAL_FRAME_FORMFACTOR_H */
