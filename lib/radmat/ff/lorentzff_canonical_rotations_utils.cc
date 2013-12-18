/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : lorentzff_canonical_rotations_utils.cc

 * Purpose :

 * Creation Date : 10-12-2013

 * Last Modified : Wed 18 Dec 2013 11:28:55 AM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "lorentzff_canonical_rotations_utils.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/utils/levi_civita.h"
#include <sstream>

namespace radmat
{

  namespace
  {
    std::string string_mom(const mom_t &p)
    {
      std::stringstream ss;
      ss << p[0] << p[1] << p[2];
      return ss.str(); 
    }

    bool same_length(const mom_t &l, const mom_t &ll)
    {
      return ( dot_mom_t(l,l) == dot_mom_t(ll,ll) ); 
    }

    double cos_theta(const mom_t &l, const mom_t &ll)
    {
      if(is_rest(l) || is_rest(ll) )
        return 0.;

      return double(dot_mom_t(l,ll))/sqrt(double(dot_mom_t(l,l) * dot_mom_t(ll,ll))); 
    }

    int volume_element(const mom_t &l, const mom_t &r)
    {
      return ( dot_mom_t(l,cross_product(l,r)) ); 
    }

  }

  // Ax = b
  bool
    check_frame_transformation(const RotationMatrix_t *A,
        const mom_t &x, 
        const mom_t &b,
        bool print_on_false)
    {
      mom_t chk = gen_mom<0,0,0>();

      for(int i = 0; i < 3; ++i)
      {
        double res(0.); 
        for(int j = 0; j < 3; ++j)
          if ( fabs( (*A)[i+1][j+1] ) > 1e-6 )
            res += (*A)[i+1][j+1] * double(x[j]); 

        chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
      }

      if( (chk[0] != b[0]) || (chk[1] != b[1]) || (chk[2] != b[2]) )
      {
        if(print_on_false)
        {
        std::cout << __func__ << ": returning false Rp = pp\n"
          << "p = " << string_mom(x) << " pp = " << string_mom(b)
          << " Rp = " << string_mom(chk)
          << " R:" << *A << std::endl; 
        }
        return false; 
      }
      return true; 
    } 

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  
  bool
    check_total_frame_transformation(const RotationMatrix_t *R,
        const mom_t &l, 
        const mom_t &r,
        const mom_t &ll,
        const mom_t &rr,
        bool print_on_false)
    {
      bool success = check_frame_transformation(R,ll,l) ;
      success &= check_frame_transformation(R,rr,r); 
      if ( print_on_false && !!! success)
      {
        std::cout << __func__ << ": rotation failure "
          << "\n l" << string_mom(l) << "  r " << string_mom(r)
          << "\nll" << string_mom(ll) << " rr " << string_mom(rr)
          << "R: " << *R << std::endl;
      }
      return success; 
    } 



  //////////////////////////////////////////////
  //////////////////////////////////////////////

  RotationMatrix_t*
    generate_frame_transformation(const mom_t &p, const mom_t &pp)
    {
      if( dot_mom_t(p,p) != dot_mom_t(pp,pp) )
      {
        std::cout << __func__ << ": error, frames are unrelated " 
          << "\np" << string_mom(p) << " pp " << string_mom(pp) 
          << std::endl;
        exit(1); 
      }

      RotationMatrix_t *R = new Tensor<double,2>( (TensorShape<2>())[4][4], 0.); 
      RotationMatrix_t *Rp = radmat::CanonicalLatticeRotationEnv::call_factory(p);
      RotationMatrix_t *Rpp = radmat::CanonicalLatticeRotationEnv::call_factory(pp);
      R->lower_index(1); 

      for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
        {
          if( fabs( (*Rp)[i][j] ) < 1e-6 )
            continue;

          for(int k = 0; k < 4; ++k)
            (*R)[i][k] += (*Rp)[i][j] * (*Rpp)[k][j]; 
        }

      delete Rp;
      delete Rpp; 

      if( !!! check_frame_transformation(R,pp,p) )
      {
        std::cout << __func__ << ": error transforming frames "
          << "\np" << string_mom(p) << " pp " << string_mom(pp) 
          << " p = Rpp\nR:" << *R << std::endl;
        delete R; 
        exit(1); 
      }

      clean_up_rot_mat(R); 

      return R; 
    }

  //////////////////////////////////////////
  //////////////////////////////////////////

  // generate an orthogonal transformation 
  RotationMatrix_t*
    generate_triad_rotation_matrix(const mom_t &l, 
        const mom_t &r,
        const mom_t &ll, 
        const mom_t &rr)
    {
      RotationMatrix_t *R = new Tensor<double,2>( (TensorShape<2>())[4][4], 0.); 
      (*R)[0][0] = 1.; 
      R->lower_index(1); 

      // if colinear then cross product is garbage
      if( colinear_momentum(l,r) && colinear_momentum(ll,rr) )
      {
        delete R; 
        return generate_frame_transformation(l,ll); 
      }

      // make some unit vectors
      itpp::Vec<double> v1,v2,w1,w2; 
      v1 = normalize(l);
      v2 = normalize(r);

      // make some unit vectors in the other frame
      w1 = normalize(ll);
      w2 = normalize(rr); 

      // use a plus minus basis to do the mapping 
      itpp::Vec<double> vp,vm,vv,vc,wp,wm,ww,wc; 
      vp = (v1 + v2)/sqrt(itpp::sum_sqr(v1 + v2));
      vm = (v1 - v2)/sqrt(itpp::sum_sqr(v1 - v2));
      vc = itpp::cross(vp,vm);

      wp = (w1 + w2)/sqrt(itpp::sum_sqr(w1 + w2));
      wm = (w1 - w2)/sqrt(itpp::sum_sqr(w1 - w2));
      wc = itpp::cross(wp,wm); 

      // outter product produces an orthogonal transformation 
      // that maps from one frame to the other (rotation matrix)
      for(int i =0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
          (*R)[i+1][j+1] += vp(i)*wp(j);
          (*R)[i+1][j+1] += vm(i)*wm(j);
          (*R)[i+1][j+1] += vc(i)*wc(j);
        }

      // should be the solution 
      //   R*l = ll;
      //   R*r = rr; 

      return R; 
    }

  //////////////////////////////////////////
  //////////////////////////////////////////


  // use zyz -- do some algebra 
  //
  //  Mathematica :
  //      Rz[t_] := RotationMatrix[t,{0,0,1}]
  //      Ry[t_] := RotationMatrix[t,{0,1,0}]
  //
  //      Eul[a_,b_,c_] := Rz[a] . Ry[b] . Rz[c]
  //
  //      Eul[a,b,c] // MatrixForm 
  //
  // NB: This is not unique 
  Hadron::CubicCanonicalRotation_t
    generate_euler_angles(const RotationMatrix_t *R)
    {
      double alpha,beta,gamma;

      if( (*R)[3][3] != 0. ) // isolates cosine of beta
      {
        beta = acos( (*R)[3][3] ); 
        if ( sin(beta) != 0. ) // can use the edge elements
        {
          double sb = sin(beta); 
          alpha = atan2( (*R)[2][3]/sb , (*R)[1][3]/sb ); 
          gamma = atan2( (*R)[3][2]/sb , -(*R)[3][1]/sb ); 
        }
        else // pick something that makes the edge elems zero and follow through
        {    // gimble lock pick a choice of gamma = pi
          gamma = acos(-1.); 
          alpha = acos( (*R)[2][2]/cos(gamma) ); 
        }
      }
      else
      { 
        // beta is pi/2 or 3pi/2 -- pick pi/2
        beta = acos(-1) / 2.;
        double sb = sin(beta); 
        alpha = atan2( (*R)[2][3]/sb , (*R)[1][3]/sb ); 
        gamma = atan2( (*R)[3][2]/sb , -(*R)[3][1]/sb ); 
      }

      Hadron::CubicCanonicalRotation_t e;
      e.alpha = alpha;
      e.beta = beta;
      e.gamma = gamma;

      return e; 
    }

  //////////////////////////////////////////////
  //////////////////////////////////////////////

  void clean_up_rot_mat(RotationMatrix_t *R)
  {
    for(int i = 1; i < 4; ++i)
      for(int j = 1; j < 4; ++j)
        if( fabs((*R)[i][j]) < 1e-6)
          (*R)[i][j] = 0.;
  }

  //////////////////////////////////////////////
  //////////////////////////////////////////////

  int dot_mom_t(const mom_t &left, const mom_t &right)
  {
    return left[0]*right[0] + left[1]*right[1] + left[2]*right[2]; 
  }

  //////////////////////////////////////////////
  //////////////////////////////////////////////

  mom_t cross_product(const mom_t &l, const mom_t &r)
  {
    mom_t c = gen_mom<0,0,0>(); 
    c[0] = l[1]*r[2] - l[2]*r[1];
    c[1] = -(l[0]*r[2] - l[2]*r[0]);
    c[2] = l[0]*r[1] - l[1]*r[0]; 
    return c; 
  } 

  //////////////////////////////////////////////
  //////////////////////////////////////////////

  bool is_rest(const mom_t &p)
  {
    return (dot_mom_t(p,p) == 0);
  }


  //////////////////////////////////////////////
  //////////////////////////////////////////////

  bool colinear_momentum(const mom_t &l, const mom_t &r)
  {
    itpp::Vec<double> nl(normalize(l)),nr(normalize(r));
    return (nl == nr || nl == -nr ); 
  }

  //////////////////////////////////////////////
  //////////////////////////////////////////////

  // test the dot product and volume element are the same
  bool related_by_rotation(const mom_t &left, const mom_t &right, 
      const mom_t &lleft, const mom_t &rright,
      const bool allow_flip)
  {
    if( is_rest(left) || is_rest(lleft) || is_rest(right) || is_rest(rright) )
    {
      std::cout << __func__ << ": passed in a rest momentum, throwing an error"
        << std::endl;
      throw std::string("related_by_rotation rest error"); 
    }

    bool success = true; 
    success &= (cos_theta(left,right) == cos_theta(lleft,rright));
    success &= (volume_element(left,right) == volume_element(lleft,rright));

    if( allow_flip ) 
    {
      if(same_length(left,lleft))
        success &= same_length(right,rright);  
      else if(same_length(left,rright)) 
        success &= same_length(lleft,right);  
      else
        return false; 
    }
    else
    {
      success &= same_length(left,lleft); 
      success &= same_length(right,rright);  
    }

    if( !!! success )
      return false; 

    // R *ll = l && R *rr = r -- doesnt work for rest
    RotationMatrix_t *Rtriad = generate_triad_rotation_matrix(left,right,lleft,rright); 
    // check that the frame transformation actually works
    success &= check_total_frame_transformation(Rtriad,left,right,lleft,rright,true); 

    // not a proper rotation 
    if( fabs( determinant(Rtriad) - 1. ) > 1e-6 )
    {
      std::cout << __func__ << ": triad had det " << determinant(Rtriad) 
        << ": which was not 1, confused Rtriad:" << *Rtriad << std::endl;
      success = false; 
    }

    delete Rtriad; 
    return success;
  }

  //////////////////////////////////////////////
  //////////////////////////////////////////////

  itpp::Vec<double> normalize(const mom_t &m)
  {
    itpp::Vec<double> r(3); 
    double norm = sqrt( double(dot_mom_t(m,m)) ); 
    r[0] = double(m[0])/norm;
    r[1] = double(m[1])/norm;
    r[2] = double(m[2])/norm;
    return r; 
  }


  // use the levi civita symbol to compute a determinant
  double determinant(const RotationMatrix_t * R)
  {
    Tensor<double,4> levi = levi_civita<double,4>(); 
    double sum; 

    for(int i = 0; i < 4; ++i)
      for(int j = 0; j < 4; ++j)
      {
        if( i == j )
          continue; 
        for(int k = 0; k < 4; ++k)
        {
          if( i == k || j == k )
            continue; 

          for(int l = 0; l < 4; ++l)
            sum += (levi[i][j][k][l] 
                * (*R)[0][i] 
                * (*R)[1][j] 
                * (*R)[2][k] 
                * (*R)[3][l]);
        } // k 
      } // j 
    return sum; 
  }

}
