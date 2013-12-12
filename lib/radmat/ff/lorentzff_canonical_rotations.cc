/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_rotations.cc

 * Purpose :


 * Last Modified : Wed 11 Dec 2013 07:57:19 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "lorentzff_canonical_rotations.h"
#include "lorentzff_canonical_rotations_utils.h"
#include "formfac/formfac_qsq.h"
#include "io/adat_xmlio.h"
#include "radmat/utils/pow2assert.h"
#include "itpp/itbase.h"
#include <sstream>
#include <exception>
#include <omp.h>
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"



namespace radmat
{


  namespace LatticeRotationEnv
  {

    namespace
    {
      void check_exit_LG(const mom_t &p, const std::string &LG)
      {
        if( Hadron::generateLittleGroup(p) != LG )
        {
          std::cout << __func__ << ": error, expected LG " 
            << LG << ", got momentum " << p[0] 
            << p[1] << p[2] << "  -->  " << Hadron::generateLittleGroup(p)
            << ": exiting" << std::endl;
          exit(1); 
        }
      }


      bool
        check_rotation_invariant(const RotationMatrix_t *R,
            const mom_t &l)
        {
          mom_t chk = gen_mom<0,0,0>();

          for(int i = 0; i < 3; ++i)
          {
            double res(0.); 
            for(int j = 0; j < 3; ++j)
              if ( fabs( (*R)[i+1][j+1] ) > 1e-6 )
                res += (*R)[i+1][j+1] * double(l[j]); 

            chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
          }

          if( (chk[0] != l[0]) || (chk[1] != l[1]) || (chk[2] != l[2]) )
          {
            std::cout << __func__ << ": Rotation error," 
              << " in " << l[0] << l[1] << l[2] 
              << " out " << chk[0] << chk[1] << chk[2]
              << "\n R:\n" << *R << std::endl;

            return false; 
          }

          return true; 
        } 

      void clean_up(RotationMatrix_t *R)
      {
        for(int i = 1; i < 4; ++i)
          for(int j = 1; j < 4; ++j)
            if( fabs((*R)[i][j]) < 1e-6)
              (*R)[i][j] = 0.;
      }

      //
      //  This is where we pick how to rotate
      //


      // the canonical frame defines the z axis  
      // everything follows as a rotaion from this frame
      RotationMatrix_t *
        rest_specialization(const mom_t &pc, const mom_t &p)
        {
          std::string LGc = Hadron::generateLittleGroup(pc);
          check_exit_LG(p,LGc);

          // if they are both at rest return identity
          if( is_rest(p) )
          {
            Tensor<double,2> *I4; 
            I4 = new Tensor<double,2>((TensorShape<2>())[4][4],0.);

            (*I4)[0][0] = 1.;
            (*I4)[1][1] = 1.;
            (*I4)[2][2] = 1.;
            (*I4)[3][3] = 1.;

            I4->lower_index(1); 

            return I4; 
          }

          RotationMatrix_t *Ref = radmat::CanonicalRotationEnv::call_factory(LGc); 
          RotationMatrix_t *Rc = radmat::CanonicalLatticeRotationEnv::call_factory(pc);
          RotationMatrix_t *Rp = radmat::CanonicalLatticeRotationEnv::call_factory(p);
          RotationMatrix_t *R= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
          RotationMatrix_t *Rtemp= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
          R->lower_index(1); 

          //
          // Rp takes me from ref to p
          // Rc takes me from ref to c 
          //
          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            {
              if( fabs( (*Rp)[i][j] ) < 1e-6 )
                continue; 

              for(int k = 0; k < 4; ++k)
                (*Rtemp)[i][k] += (*Rp)[i][j] * (*Rc)[k][j]; 
            }

          // hit it with a ref rotation  
          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            {
              if( fabs( (*Rtemp)[i][j] ) < 1e-6 )
                continue; 

              for(int k = 0; k < 4; ++k)
                (*R)[i][k] += (*Rtemp)[i][j] * (*Ref)[k][j]; 
            }

          delete Ref; 
          delete Rc; 
          delete Rp; 
          delete Rtemp;

          return R; 
        }


      RotationMatrix_t *
        generate_right_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          std::string LGl = Hadron::generateLittleGroup(cl);
          std::string LGr = Hadron::generateLittleGroup(cr);

          check_exit_LG(l,LGl);
          check_exit_LG(r,LGr); 

          // if righty is at rest then he needs 
          // to follow lefty
          if( is_rest(cr) ) 
            return rest_specialization(cl,l); 

          // righty doesn't care if lefty is at rest, that
          // is lefty's problem 

          // return the normal rotations with the ref built in
          return radmat::CanonicalRotationEnv::call_factory(r);
        }



      RotationMatrix_t *
        generate_left_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          std::string LGl = Hadron::generateLittleGroup(cl);
          std::string LGr = Hadron::generateLittleGroup(cr);

          check_exit_LG(l,LGl);
          check_exit_LG(r,LGr); 

          // if righty is at rest then righty must 
          // follow lefty 
          if( is_rest(cr) )
            return radmat::CanonicalRotationEnv::call_factory(l); 

          // if lefty is at rest then we need to follow righty 
          if( is_rest(cl) )
            return rest_specialization(cr,r); 

          // need to be deleted
          // 
          //    These are lattice rotations 
          //
          RotationMatrix_t *Rcl = radmat::CanonicalLatticeRotationEnv::call_factory(cl);
          RotationMatrix_t *Rcr = radmat::CanonicalLatticeRotationEnv::call_factory(cr);
          RotationMatrix_t *Rl = radmat::CanonicalLatticeRotationEnv::call_factory(l);
          RotationMatrix_t *Rr = radmat::CanonicalLatticeRotationEnv::call_factory(r);
          RotationMatrix_t *Rtemp= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
          RotationMatrix_t *Rcr_r= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
          RotationMatrix_t *Rcl_l= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
          Rtemp->lower_index(1); 

          // gets returned
          RotationMatrix_t *R= new Tensor<double,2>((TensorShape<2>())[4][4],0.);
          R->lower_index(1);

          clean_up(Rcl);
          clean_up(Rcr); 
          clean_up(Rl);
          clean_up(Rr); 

          // Rcr_r is the rotation that takes me from the 
          // right reference frame to the right frame 
          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            {
              if( fabs( (*Rr)[i][j] ) < 1e-6 )
                continue; 

              for(int k = 0; k < 4; ++k)
                (*Rcr_r)[i][k] += (*Rr)[i][j] * (*Rcr)[k][j];   
            }


          // Rcl_l is the rotation that takes me from the 
          // left reference frame to the left frame 
          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            {
              if( fabs( (*Rl)[i][j] ) < 1e-6 )
                continue; 

              for(int k = 0; k < 4; ++k)
                (*Rcl_l)[i][k] += (*Rl)[i][j] * (*Rcl)[k][j];   
            }

          clean_up(Rcr_r); 
          clean_up(Rcl_l); 

          bool frame_tform = true; 
          frame_tform &= check_frame_transformation(Rcr_r,cr,r);
          frame_tform &= check_frame_transformation(Rcl_l,cl,l); 

          if( !!! frame_tform )
          {
            std::cout << __func__ << ": frame_transformatio error"
              << std::endl;
            exit(1); 
          }

          // Rtemp is is the rotation that takes me from the 
          //  left frame to the left reference frame 
          //  followed by the rotation that takes me 
          //  from the right reference frame to the right frame 
          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
            {
              if( fabs( (*Rcr_r)[i][j] ) < 1e-6 )
                continue;

              for(int k = 0; k < 4; ++k)
                (*Rtemp)[i][k] += (*Rcr_r)[i][j] * (*Rcl_l)[k][j];
            }


          // this can only rotate about the momentum
          if( !!! check_rotation_invariant(Rtemp,l))
          {
            std::cout << "\n\n" << __func__ << ": cl" << cl[0] << cl[1] << cl[2] 
              << "  l" << l[0] << l[1] << l[2] 
              << "  cr" << cr[0] << cr[1] << cr[2] 
              << "  r" << r[0] << r[1] << r[2] 
              << "\n*************************"
              << std::endl;

            clean_up(Rcl); 
            clean_up(Rcr);
            clean_up(Rl);
            clean_up(Rr);

            std::cout << "\nRcl:" << *Rcl << "\nRcr:" << *Rcr
              << "\nRl:" << *Rl << "\nRr:" << *Rr 
              << "\nRcr_r:" << *Rcr_r << "\nRcl_l:" << *Rcl_l 
              << std::endl;

            exit(1); 
          } 

          delete Rl;

          // stick the left rotation (and ref) back in 
          Rl = radmat::CanonicalRotationEnv::call_factory(l); 

          // put the Rl back in 
          for(int i = 0; i < 4; ++i)
            for(int j = 0; j < 4; ++j)
              for(int k = 0; k < 4; ++k)
                (*R)[i][k] += (*Rtemp)[i][j] * (*Rl)[j][k]; 

          clean_up(R); 

          delete Rtemp;
          delete Rcl;
          delete Rcr;
          delete Rl;
          delete Rr; 

          return R; 
        }


    } // anonomyous




    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////

    namespace
    {
      bool local_registration = false; 

      rHandle<RotationMatrix_t>
        left_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          return rHandle<RotationMatrix_t>(generate_left_rotation(cl,cr,l,r)); 
        }

      rHandle<RotationMatrix_t>
        right_rotation(const mom_t &cl, 
            const mom_t &cr, 
            const mom_t &l, 
            const mom_t &r) 
        {
          return rHandle<RotationMatrix_t>(generate_right_rotation(cl,cr,l,r)); 
        }

    } // anonomyous


    bool registerAll(void)
    {
      bool success = true; 

      if ( !!! local_registration)
      {
        success &= TheRotationGroupGenerator::Instance().BlastOff(); 
        local_registration = true; 
      }

      if( !!! success )
      {
        throw std::string("reg error in LatticeRotationEnv");
      }

      return success; 
    }

    std::string
      rotation_group_label(const mom_t &l, const mom_t &r)
      {
        return TheRotationGroupGenerator::Instance().get_can_frame_string(l,r);  
      }

    rHandle<RotationMatrix_t> 
      get_left_rotation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can; 
        can = TheRotationGroupGenerator::Instance().get_can_frame(l,r); 
        return left_rotation(can.first,can.second,l,r); 
      }

    rHandle<RotationMatrix_t> 
      get_right_rotation(const mom_t &l, const mom_t &r)
      {
        std::pair<mom_t,mom_t> can; 
        can = TheRotationGroupGenerator::Instance().get_can_frame(l,r); 
        return right_rotation(can.first,can.second,l,r); 
      }

  } // LatticeRotationEnv

} // radmat






//      // this is exact, it basically cooks up an orthogonal transformation
//      //
//      //    -- we will get screwed and have to play a phase game 
//      itpp::Mat<double> 
//        triad_rotation_matrix(const mom_t &c_left, 
//            const mom_t &c_right,
//            const mom_t &left, 
//            const mom_t &right)
//        {
//          itpp::Mat<double> A(3,3),B(3,3); 
//
//          A.set_row(0,normalize(c_left)); 
//          A.set_row(1,normalize(c_right)); 
//          A.set_row(3,normalize(cross_product(c_left,c_right))); 
//
//          B.set_row(0,normalize(left)); 
//          B.set_row(1,normalize(right)); 
//          B.set_row(2,normalize(cross_product(left,right))); 
//
//          return A * itpp::transpose(B); 
//        }
//
//
//      template<int A, int B, int C, int X, int Y, int Z> 
//        RotationMatrix_t * 
//        generate_rotation(const mom_t & left, const mom_t &right)
//        {
//          RotationMatrix_t * R;
//          R = new RotationMatrix_t( (TensorShape<2>())[4][4] , 0. );  
//
//          mom_t cleft = gen_mom<A,B,C>(); 
//          mom_t cright = gen_mom<X,Y,Z>();  
//
//          // sanity
//          POW2_ASSERT( related_by_rotation(left,right,cleft,cright) );          
//
//          itpp::Mat<double> R3 = triad_rotation_matrix(cleft,cright,left,right); 
//
//          (*R)[0][0] = 1.;
//
//          for(int i = 0; i < 3; ++i)
//            for(int j = 0; j < 3; ++j)
//              (*R)[i+1][j+1] = R3(i,j); 
//
//          return R; 
//        }

