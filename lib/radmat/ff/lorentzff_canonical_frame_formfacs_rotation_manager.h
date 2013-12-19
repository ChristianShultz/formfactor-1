#ifndef LORENTZFF_CANONICAL_FRAME_FORMFACS_ROTATION_MANAGER_H
#define LORENTZFF_CANONICAL_FRAME_FORMFACS_ROTATION_MANAGER_H 

#include "lorentzff_Wigner_D_matrix_factory.h"
#include "lorentzff_canonical_rotations_utils.h"
#include "lorentzff_canonical_rotations.h"
#include "lorentzff_formfac_utils.h"
#include "ff_gen_llsq_row.h"
#include "hadron/irrep_util.h"
#include "radmat/utils/tensor.h"
#include <complex>
#include <sstream>

namespace radmat
{

  // most positive helicity is mapped to zero
  struct DMatrixManager
  {
    typedef radmat::LatticeRotationEnv::FrameOrientation_t
      FrameOrientation_t; 

    virtual ~DMatrixManager(void) {}

    virtual std::string 
      string_mom(const mom_t &p) const
    {
      std::stringstream ss; 
      ss << p[0] << p[1] << p[2] ;
      return ss.str(); 
    }

    virtual FrameOrientation_t 
      get_frame(const mom_t &l, 
          const mom_t &r) const ;

    virtual void 
      check_throw_frame_err(const RotationMatrix_t* Rtriad, 
          const FrameOrientation_t &) const; 

    virtual WignerMatrix_t* 
      get_can_mat(const mom_t &p, 
          const int J) const;

    virtual void 
      conjugate(WignerMatrix_t * D) const;

    virtual void 
      transpose(WignerMatrix_t *D) const; 

    virtual void
      dagger(WignerMatrix_t *D) const; 

    virtual void 
      clean(WignerMatrix_t *D, 
          const double thresh=1e-6) const;  

    virtual RotationMatrix_t*
      triad_rotation_matrix(const mom_t &l, 
          const mom_t &r) const; 

    virtual WignerMatrix_t*
      triad_rotation_wigner_matrix(const RotationMatrix_t *R, 
          const mom_t &l, 
          const mom_t &r, 
          const int J) const;

    virtual WignerMatrix_t*
      left_wigner_matrix(const RotationMatrix_t *R, 
          const mom_t &l, 
          const mom_t &r, 
          const int J,
          bool print=false) const; 

    virtual WignerMatrix_t*
      right_wigner_matrix(const RotationMatrix_t *R,
          const mom_t &l, 
          const mom_t &r, 
          const int J,
          bool print = false) const; 
  }; 



  template<class DerivedFF>
    struct FormFacRotationManager
    : public DMatrixManager
  {
    typedef DMatrixManager::FrameOrientation_t FrameOrientation_t; 
    typedef Tensor<double,1> p4_t; 

    virtual ~FormFacRotationManager() {}

    virtual std::pair<mom_t,mom_t> 
      pair_mom(const p4_t &l, const p4_t &r, const double kick) const
      {
        return std::pair<mom_t,mom_t>(get_space_mom(l,kick),get_space_mom(r,kick)); 
      }


    virtual std::pair<p4_t,p4_t>
      can_mom(const RotationMatrix_t *R,
          const p4_t &l, 
          const p4_t &r,
          const double kick) const
      {
        Tensor<double,1> ll( (TensorShape<1>())[4], 0.);
        Tensor<double,1> rr( (TensorShape<1>())[4], 0.);

        // rotate by the transpose
        for( int i = 0; i < 4; ++i)
          for(int j = 0; j < 4; ++j)
          {
            ll[j] += (*R)[i][j] * l[i]; 
            rr[j] += (*R)[i][j] * r[i];
          }

        std::pair<mom_t,mom_t> chk = pair_mom(ll,rr,kick); 
        std::pair<mom_t,mom_t> fmom = pair_mom(l,r,kick); 

        // are these the momentum i think that they are? 
        if ( !!! check_total_frame_transformation( R, fmom.first,fmom.second,chk.first,chk.second,true ) )
        {
          std::cout << __func__ << ": frame transformation error" << std::endl;
          
          std::cout << "int moms l" << string_mom(fmom.first) << " r " << string_mom(fmom.second)
            << " ll " << string_mom(chk.first) << " rr " << string_mom(chk.second) << std::endl;

          std::cout << "pl " << l << " pll " << ll << " pr " << r << " prr "
            << rr << " R " << *R << std::endl;
          throw std::string("frame remap error"); 
        }

        return std::pair<p4_t,p4_t>(ll,rr); 
      }

    virtual Tensor<std::complex<double> , 1> 
      rotate(const RotationMatrix_t *R, 
          const Tensor<std::complex<double>,1> &in) const
      {
        Tensor<std::complex<double>,1> out( (TensorShape<1>())[4] , std::complex<double>(0.,0.) ); 

        for(int i = 0; i < 4; ++i)
          for(int j = 0; j < 4; ++j)
            out[i] += (*R)[i][j] * in[j];

        return out; 
      }

    virtual Tensor<std::complex<double>,1>
      operator()(const p4_t &l, 
          const p4_t &r, 
          const double kick, 
          const int Jl, 
          const int Jr, 
          const int hl, 
          const int hr) const
      {
        std::pair<mom_t,mom_t> moms = pair_mom(l,r,kick);  
        RotationMatrix_t *R = triad_rotation_matrix(moms.first,moms.second); 
        std::pair<p4_t,p4_t> can_moms = can_mom(R,l,r,kick); 
        WignerMatrix_t * Wl = left_wigner_matrix(R,moms.first,moms.second,Jl);
        WignerMatrix_t * Wr = right_wigner_matrix(R,moms.first,moms.second,Jr);

        Tensor<std::complex<double>,1> sum( (TensorShape<1>())[4] , std::complex<double>(0.,0.) ); 

        int left_h = Jl - hl; 
        int right_h = Jr - hr; 
        int left_bound = 2*Jl +1; 
        int right_bound = 2*Jr +1; 

        for(int lh = 0; lh < left_bound; ++lh)
          for(int rh = 0; rh < right_bound; ++rh)
            sum += (
                ( (*Wl)[left_h][lh] * (*Wr)[rh][right_h] ) // D matrix phase
                * impl(can_moms.first,can_moms.second,kick,Jl-lh,Jr-rh) ); // impl in the canonical frame

        Tensor<std::complex<double>,1> ret = rotate(R,sum); 

        delete Wl;
        delete Wr; 
        delete R; 

        return ret; 
      }

    // derived classes implement impl 
    virtual Tensor<std::complex<double>, 1> 
      impl(const p4_t &l, const p4_t &r, const double kick, const int hl, const int hr) const
      {
        return static_cast<const DerivedFF*>(this)->impl(l,r,kick,hl,hr); 
      } 

    virtual std::string 
      ff_impl(void) const 
      {
        return static_cast<const DerivedFF*>(this)->ff_impl(); 
      }

  }; 


} // radmat



#endif /* LORENTZFF_CANONICAL_FRAME_FORMFACS_ROTATION_MANAGER_H */