/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_frame_formfacs_rotation_manager.cc

 * Purpose :

 * Creation Date : 18-12-2013

 * Last Modified : Wed 18 Dec 2013 02:13:49 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "lorentzff_formfac_utils.h"
#include "hadron/clebsch.h"



namespace radmat
{

  DMatrixManager::FrameOrientation_t 
    DMatrixManager::get_frame(const mom_t &l, const mom_t &r)
    {
      return radmat::LatticeRotationEnv::get_frame_orientation(l,r); 
    }

      // there is a delta function that MUST be satisfied 
  void 
    DMatrixManager::check_throw_frame_err(const RotationMatrix_t *R, const FrameOrientation_t &f) const
    {
      if( !!! check_total_frame_transformation(R,f.l,f.r,f.ll,f.rr,true) )
      {
        std::cout << __func__ << ": throwing string " << std::endl;
        throw std::string("triad wigner rotation error"); 
      }
    }

  WignerMatrix_t* 
    DMatrixManager::get_can_mat(const mom_t &p, const int J)
    {
      return radmat::WignerDMatrixEnv::call_factor(p,J); 
    }  

  void 
    DMatrixManager::conjugate(WignerMatrix_t * D)
    {
      WignerMatrix_t::iterator it;
      for(it = D->begin(); it != D->end(); ++it)
        *it = std::conj(*it); 
    }

  void 
    DMatrixManager::clean(WignerMatrix_t *D, const double thresh)
    {
      WignerMatrix_t::iterator it; 
      for(it = D->begin(); it != D->end(); ++it)
        *it = round_to_zero( *it , thresh ); 
    }

  RotationMatrix_t*
    DMatrixManager::triad_rotation_matrix(const mom_t &l, const mom_t &r) const
    {
      FrameOrientation_t f = get_frame(l,r); 
      return generate_triad_rotation_matrix(f.l,f.r,f.ll,f.rr); 
    }

  WignerMatrix_t*
    DMatrixManager::triad_rotation_wigner_matrix(const RotationMatrix_t *R, const mom_t &l, const mom_t &r, const int J)
    {
      FrameOrientation_t f = get_frame(l,r); 
      check_throw_frame_err(R,f); 

      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   

      Hadron::CubicCanonicalRotation_t eul = generate_euler_angles(R); 

      for(int m1 = -J; m1 <= J; ++m1)
        for(int m2 = -J; m2 <= J; ++m2)
        {
          std::complex<double> cd = SEMBLE::toScalar(
              Hadron::Wigner_D(2*J,2*m1,2*m2,eul.alpha,eul.beta,eul.gamma));
          (*W)[J-m1][J-m2] = round_to_zero(cd,1e-6); 
        }

      return W; 
    }


  // notation follows notes
  WignerMatrix_t* 
    DMatrixManager::left_wigner_matrix(const RotationMatrix_t *R, const mom_t &l, const mom_t &r, const int J) const
    {
      FrameOrientation_t f = get_frame(l,r); 
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   
      WignerMatrix_t *Wt,*Wn,*Wi; 

      // the delta function is checked here
      Wt = triad_rotation_wigner_matrix(R,l,r,J); 
      Wn = radmat::WignerDMatrixEnv::call_factory(f.l,J);
      Wi = radmat::WignerDMatrixEnv::call_factory(f.cl,J);

      conjugate(Wn); 

      for(int i = 0; i < bound; ++i)
        for(int j = 0; j < bound; ++j)
          for(int k = 0; k < bound; ++k)
            for(int l = 0; l < bound; ++l)
              (*W)[i][l] += (*Wn)[i][j] * (*Wt)[j][k] * (*Wi)[k][l];

      clean(W); 

      delete Wt;
      delete Wn;
      delete Wi; 

      return W; 
    }

  // notation follows notes
  WignerMatrix_t*
    DMatrixManager::right_wigner_matrix(const RotationMatrix_t *R, const mom_t &l, const mom_t &r, const int J) const
    {
      FrameOrientation_t f = get_frame(l,r); 
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   
      WignerMatrix_t *Wt,*Wk,*Wl; 

      // the delta function is checked here
      Wt = triad_rotation_wigner_matrix(R,l,r,J); 
      Wl = radmat::WignerDMatrixEnv::call_factory(f.r,J);
      Wk = radmat::WignerDMatrixEnv::call_factory(f.cr,J);

      conjugate(Wk);
      conjugate(Wt); 

      for(int i = 0; i < bound; ++i)
        for(int j = 0; j < bound; ++j)
          for(int k = 0; k < bound; ++k)
            for(int l = 0; l < bound; ++l)
              (*W)[i][l] += (*Wk)[i][j] * (*Wt)[j][k] * (*Wl)[k][l];

      clean(W); 

      delete Wt;
      delete Wl;
      delete Wk; 

      return W; 
    }



}// radmat
