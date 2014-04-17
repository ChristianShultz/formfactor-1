/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : Wigner_D_matrix_manager.cc

 * Purpose :

 * Creation Date : 14-04-2014

 * Last Modified : Thu 17 Apr 2014 11:36:58 AM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "Wigner_D_matrix_manager.h"
#include "hadron/clebsch.h"
#include "radmat/utils/pow2assert.h"
#include "semble/semble_semble.h"
#include "rotation_group_generator.h"
#include "radmat/utils/printer.h"
#include "Wigner_D_matrix_factory.h"


namespace radmat
{
  namespace
  {
    std::complex<double> complex_zero(0.,0.); 

    std::complex<double> round_to_zero(const std::complex<double> &cd, const double thresh=1e-6)
    {return ( std::norm(cd) < thresh ) ? complex_zero : cd ; }


    // a momentum key 
    struct wig_key
    {
      wig_key() {}
      wig_key(const int xx, const int yy, const int zz)
        : x(xx) , y(yy) , z(zz)
      { }

      wig_key(const ADATXML::Array<int> &p)
        : x(p[0]) , y(p[1]), z(p[2]) 
      { }

      int x,y,z; 
    }; 

    // error streaming 
    std::ostream & operator<<(std::ostream &o, const wig_key &k)
    {return (o << "p" << k.x << k.y << k.z);}

    // a frame key 
    struct wig_pair_key
    {
      wig_pair_key() {}
      wig_pair_key(const ADATXML::Array<int> &ll , 
          const ADATXML::Array<int> &rr, 
          const int JJ)
        : l(ll) , r(rr), J(JJ)
      {}

      wig_key l,r;
      int J;  
    }; 

    // error streaming 
    std::ostream &operator<<(std::ostream &o, const wig_pair_key &k)
    { return (o << "l" << k.l << "r" << k.r << "J" << k.J);}


    // key comparison class, bung them all together  
    struct wig_key_comp
    {
      bool operator()(const wig_key &l, const wig_key &r) const
      {
        if(l.x != r.x)
          return l.x < r.x; 
        if(l.y != r.y)
          return l.y < r.y; 
        return l.z < r.z; 
      }

      bool operator()(const wig_pair_key &l, const wig_pair_key &r) const
      {
        if ( !!! exact_equivalence(l.l,r.l) )
          return this->operator()(l.l,r.l);

        if( !!! exact_equivalence(l.r,r.r) )
          return this->operator()(l.r,r.r); 

        return l.J < r.J; 
      }

      bool exact_equivalence( const wig_key &l, const wig_key &r) const 
      {
        return ((l.x==r.x) && (l.y==r.y) && (l.z==r.z));
      }

      bool exact_equivalence(const wig_pair_key &l , const wig_pair_key &r) const
      {
        return exact_equivalence(l.l,r.l) && exact_equivalence(l.r,r.r) && (l.J == r.J); 
      }

    };


    typedef std::map<wig_pair_key,WignerMatrix_t,wig_key_comp> WignerMatrixMap_t; 
    WignerMatrixMap_t left_map, right_map; 

    WignerMatrixMap_t* call_left_map()
    {
      return &left_map;
    }

    WignerMatrixMap_t* call_right_map()
    {
      return &right_map;
    }

    struct injection_map_printer
    {
      static void print(const std::string &msg)
      {}
      // {std::cout << "injection_map_printer " << msg << std::endl;}
    };

    // pre compute avaliable wigner matricies  
    void inject_map(const int J)
    {
      typedef radmat::LatticeRotationEnv::TheRotationGroupGenerator RG; 
      std::map<mom_pair_key,mom_pair_key,mom_key_comp>::const_iterator it; 

      DMatrixManager Wigner; 

      for(it = RG::Instance().can_frame_map.begin(); it != RG::Instance().can_frame_map.end(); ++it)
      {
        mom_t l, r; 
        l = it->first.l.mom(); 
        r = it->first.r.mom(); 
        RotationMatrix_t * R = Wigner.rotation_matrix(l,r); 
        WignerMatrix_t * left = Wigner.left_wigner_matrix(R,l,r,J); 
        WignerMatrix_t *right = Wigner.right_wigner_matrix(R,l,r,J); 
        wig_pair_key key(l,r,J); 

        //  std::stringstream ss; 
        //  ss << key; 
        //  printer_function<injection_map_printer>(ss.str()); 

        call_left_map()->insert(std::make_pair(key,*left)); 
        call_right_map()->insert(std::make_pair(key,*right)); 

        delete R; 
        delete left; 
        delete right; 
      }
    }




    bool local_registration = false;  
    bool use_wigner_map = false; 

    bool init_maps(const int J)
    {
      if( !!! local_registration )
      {
        for(int i = 0; i <= J; ++i)
          inject_map(i); 

        local_registration = true; 
        use_wigner_map = true; 
      }
      return true; 
    }

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
        bool use_map) const
    {
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   

      // use the precomputed wigner matricies
      if ( use_wigner_map && use_map)
      {
        wig_pair_key key(l,r,J); 
        WignerMatrixMap_t::iterator it = call_left_map()->find(key); 
        if( it == call_left_map()->end() )
        {
          std::cout << __PRETTY_FUNCTION__ << ": error missing key " 
            << key << std::endl;
          exit(1); 
        }

        *W = it->second; 
      }
      else
      {
        std::pair<mom_t,mom_t> can = get_frame(l,r); 
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


        delete Wt;
        delete Wn;
        delete Wi; 
      }

      return W; 
    }

  // notation follows notes
  WignerMatrix_t*
    DMatrixManager::right_wigner_matrix(const RotationMatrix_t *R, 
        const mom_t &l, 
        const mom_t &r, 
        const int J,
        bool use_map) const
    {
      int bound = 2*J+1; 
      WignerMatrix_t *W  = new WignerMatrix_t( (TensorShape<2>())[bound][bound] , std::complex<double>(0.,0.));   

      // use the precomputed wigner matricies
      if( use_wigner_map && use_map )
      {
        wig_pair_key key(l,r,J); 
        WignerMatrixMap_t::iterator it = call_right_map()->find(key); 

        if( it == call_left_map()->end() )
        {
          std::cout << __PRETTY_FUNCTION__ << ": error missing key " 
            << key << std::endl;
          exit(1); 
        }

        *W =  it->second; 

      }
      else
      {
        std::pair<mom_t,mom_t> can = get_frame(l,r); 
        int bound = 2*J+1; 
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

        delete Wt;
        delete Wl;
        delete Wk; 

      }

      return W; 
    }


  namespace WignerThreadMapEnv
  {
    bool registerAll(){ return init_maps(4); }
  }



}// radmat

