/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : Wigner_D_matrix_manager.cc

* Purpose :

* Creation Date : 14-04-2014

* Last Modified : Mon 14 Apr 2014 05:43:24 PM EDT

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "Wigner_D_matrix_manager.h"
#include "hadron/clebsch.h"
#include "radmat/utils/pow2assert.h"
#include "semble/semble_semble.h"
#include "rotation_group_generator.h"



namespace radmat
{
  namespace
  {
    std::complex<double> complex_zero(0.,0.); 

    std::complex<double> round_to_zero(const std::complex<double> &cd, const double thresh=1e-6)
    {return ( std::norm(cd) < thresh ) ? complex_zero : cd ; }
  } // anonomyous 
    

  std::pair<mom_t,mom_t>
    DMatrixManager::get_frame(const mom_t &l, const mom_t &r) const
    {
      return radmat::LatticeRotationEnv::rotation_group_key(l,r); 
    }

      // there is a delta function that MUST be satisfied 
  void 
    DMatrixManager::check_throw_frame_err(const RotationMatrix_t *R, 
        const std::pair<mom_t,mom_t> &f, 
        const std::pair<mom_t,mom_t> &c) const
    {
      if( !!! check_total_frame_transformation(R,f.first,f.second,c.first,c.second,true) )
      {
        std::cout << __func__ << ": throwing string " << std::endl;
        throw std::string("triad wigner rotation error"); 
      }
    }

  WignerMatrix_t* 
    DMatrixManager::get_can_mat(const mom_t &p, const int J) const
    {
      return radmat::WignerDMatrixEnv::call_factory(p,J); 
    }  

  void 
    DMatrixManager::conjugate(WignerMatrix_t * D) const
    {
      WignerMatrix_t::iterator it;
      for(it = D->begin(); it != D->end(); ++it)
        *it = std::conj(*it); 
    }

  void 
    DMatrixManager::transpose(WignerMatrix_t *D) const
    {
      WignerMatrix_t foo(*D); 
      std::vector<idx_t> dimensions = D->getDim(); 
      POW2_ASSERT( dimensions.size() == 2 );
      POW2_ASSERT( dimensions[0] == dimensions[1] );
      int bound = dimensions[0]; 
      for(int i = 0; i < bound; ++i)
        for(int j =0; j < bound; ++j)
          (*D)[i][j] = foo[j][i]; 
    }

  void 
    DMatrixManager::dagger(WignerMatrix_t *D) const
    {
      conjugate( D ); 
      transpose( D ); 
    }


  void 
    DMatrixManager::clean(WignerMatrix_t *D, const double thresh) const
    {
      WignerMatrix_t::iterator it; 
      for(it = D->begin(); it != D->end(); ++it)
        *it = round_to_zero( *it , thresh ); 
    }

  RotationMatrix_t*
    DMatrixManager::rotation_matrix(const mom_t &l, const mom_t &r) const
    {
      std::pair<mom_t,mom_t> f  = get_frame(l,r); 
      return generate_rotation_matrix(l,r,f.first,f.second); 
    }

  WignerMatrix_t*
    DMatrixManager::wigner_matrix(const RotationMatrix_t *R,
        const mom_t &l,
        const mom_t &r,
        const int J) const
    {
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
    DMatrixManager::left_wigner_matrix(const RotationMatrix_t *R,
        const mom_t &l,
        const mom_t &r, 
        const int J,
        bool print) const
    {
      std::pair<mom_t,mom_t> can = get_frame(l,r); 
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   
      WignerMatrix_t *Wt,*Wn,*Wi; 

      // the delta function is checked here
      Wt = wigner_matrix(R,l,r,J); 
      Wn = radmat::WignerDMatrixEnv::call_factory(l,J);
      Wi = radmat::WignerDMatrixEnv::call_factory(can.first,J);

      dagger(Wi); 
      dagger(Wt); 

      for(int i = 0; i < bound; ++i)
        for(int j = 0; j < bound; ++j)
          for(int k = 0; k < bound; ++k)
            for(int l = 0; l < bound; ++l)
              (*W)[i][l] += (*Wi)[i][j] * (*Wt)[j][k] * (*Wn)[k][l];

      dagger(W); 
      clean(W); 

      if( print ) 
      {
        std::cout << __func__ << ": moms " << "l" << string_mom(l)
          << " r " << string_mom(r) << " cl " << string_mom(can.first)
          << " cr " << string_mom(can.second) << std::endl;

        clean(Wn); 
        clean(Wt);
        clean(Wi); 

        std::cout << __func__ << ": ingredients "  
          << "Wn:" <<  *Wn << "\nWt:" << *Wt << "\nWi:" << *Wi
          << std::endl;
      }



      delete Wt;
      delete Wn;
      delete Wi; 

      return W; 
    }

  // notation follows notes
  WignerMatrix_t*
    DMatrixManager::right_wigner_matrix(const RotationMatrix_t *R, 
        const mom_t &l, 
        const mom_t &r, 
        const int J,
        bool print) const
    {
      std::pair<mom_t,mom_t> can = get_frame(l,r); 
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   
      WignerMatrix_t *Wt,*Wk,*Wl; 

      // the delta function is checked here
      Wt = wigner_matrix(R,l,r,J); 
      Wl = radmat::WignerDMatrixEnv::call_factory(r,J);
      Wk = radmat::WignerDMatrixEnv::call_factory(can.second,J);

      dagger(Wl); 

      for(int i = 0; i < bound; ++i)
        for(int j = 0; j < bound; ++j)
          for(int k = 0; k < bound; ++k)
            for(int l = 0; l < bound; ++l)
              (*W)[i][l] += (*Wl)[i][j] * (*Wt)[j][k] * (*Wk)[k][l];

      dagger(W); 
      clean(W); 


      if( print ) 
      {
        std::cout << __func__ << ": moms " << "l" << string_mom(l)
          << " r " << string_mom(r) << " cl " << string_mom(can.first)
          << " cr " << string_mom(can.second) << std::endl;

        clean(Wk); 
        clean(Wt);
        clean(Wl); 

        std::cout << __func__ << ": ingredients "  
          << "Wk:" <<  *Wk << "\nWt:" << *Wt << "\nWl:" << *Wl
          << std::endl;
      }


      delete Wt;
      delete Wl;
      delete Wk; 

      return W; 
    }



}// radmat

