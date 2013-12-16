#ifndef LORENTZFF_WIGNER_D_MATRIX_EMBEDDING_H
#define LORENTZFF_WIGNER_D_MATRIX_EMBEDDING_H 

#include "lorentzff_Wigner_D_matrix_factory.h"
#include "lorentzff_canonical_rotations.h"
#include "lorentzff_formfac_utils.h"

namespace radmat
{

  struct primitiveEmbededWignerDMatrix
  {
    virtual ~primitiveEmbededWignerDMatrix() {}

    virtual radmat::LatticeRotationEnv::FrameOrientation_t 
      get_frame(const mom_t &l, const mom_t &r) const
      {
        return  radmat::LatticeRotationEnv::get_frame_orientation(l,r); 
      }

    virtual void 
      conjugate(WignerMatrix_t* W) const
      {
        WignerMatrix_t::iterator it; 
        for(it = W->begin(); it != W->end(); ++it)
          *it = std::conj(*it); 
      }

    virtual  WignerMatrix_t
      composite_wigner_matrix(const mom_t &l, const mom_t &r) const = 0; 

    virtual void 
      clean_up(WignerMatrix_t &W) const
      {
        WignerMatrix_t::iterator it; 
        for(it = W.begin(); it != W.end(); ++it)
          *it = round_to_zero(*it,1e-6); 
      }
  };

  //////////////////////////////////////////////////////

  template<int J>
    struct embededWignerDMatrix
    :  public primitiveEmbededWignerDMatrix 
    {

      virtual ~embededWignerDMatrix() {}

      virtual  WignerMatrix_t
        composite_wigner_matrix(const mom_t &l, const mom_t &r) const
        {
          radmat::LatticeRotationEnv::FrameOrientation_t frame;
          frame = get_frame(l,r); 

          WignerMatrix_t *Dcl,*Dcr,*Dl,*Dr;
          WignerMatrix_t D( (TensorShape<2>())[2*J+1][2*J+1], std::complex<double>(0.,0.) ); 

          Dcl = radmat::WignerDMatrixEnv::call_factory<J>(frame.cl);
          Dcr = radmat::WignerDMatrixEnv::call_factory<J>(frame.cr);
          Dl = radmat::WignerDMatrixEnv::call_factory<J>(frame.l);
          Dr = radmat::WignerDMatrixEnv::call_factory<J>(frame.r);

          conjugate(Dcl);
          conjugate(Dr); 

          for(int i = 0; i < 2*J+1; ++i)
            for(int j = 0; j < 2*J+1; ++j)
              for(int k = 0; k < 2*J+1; ++k)
                for(int l = 0; l < 2*J+1; ++l)
                  for(int m = 0; m < 2*J+1; ++m)
                    D[i][m] += ( (*Dcl)[j][i] * (*Dcr)[j][k]
                        * (*Dr)[k][l] * (*Dl)[l][m] );

          delete Dcl;
          delete Dcr;
          delete Dl;
          delete Dr; 

          clean_up(D); 

          return D; 
        } 

    }; 


} // radmat
#endif /* LORENTZFF_WIGNER_D_MATRIX_EMBEDDING_H */
