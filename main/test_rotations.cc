/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : test_rotations.cc

* Purpose :

* Creation Date : 11-12-2013

* Last Modified : Thu 19 Dec 2013 09:57:43 AM EST

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/

#include "radmat/register_all/register_all.h"
#include <map>
#include <string>
#include <sstream>
#include "radmat/redstar_interface/redstar_canonical_rotations.h"
#include "radmat/redstar_interface/redstar_canonical_lattice_rotations.h"
#include "radmat/ff/lorentzff_canonical_rotations_utils.h"
#include "radmat/ff/lorentzff_canonical_frame_formfacs_rotation_manager.h"
#include "radmat/ff/lorentzff_canonical_rotations.h"
#include "radmat/ff/lorentzff_Wigner_D_matrix_factory.h"
#include "radmat/ff/lorentzff_Wigner_D_matrix_embedding.h"
#include "radmat/utils/levi_civita.h"
#include "hadron/irrep_util.h"
#include "formfac/formfac_qsq.h"

typedef radmat::mom_t mom_t; 
typedef radmat::RotationMatrix_t RotationMatrix_t; 

// a hacky thingy
struct KeyRetriever 
{
  template<typename T>
    typename T::first_type operator() (T KeyValPair) const
    {
      return KeyValPair.first; 
    }
};

/////////////////////////////////////////////
std::string string_mom(const mom_t &p)
{
  std::stringstream ss; 
  ss << p[0] << p[1] << p[2];
  return ss.str(); 
}

/////////////////////////////////////////////
mom_t ref_mom(const mom_t &p)
{
  mom_t foo = FF::canonicalOrder(p); 

  if( (foo[0] == 1) && (foo[1] == 1) && (foo[2] == 1) )
  {
    return radmat::gen_mom<1,1,1>();
  }
  else if( (foo[0] == 1) && (foo[1] == 1) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,1,1>();
  }
  else if( (foo[0] == 1) && (foo[1] == 0) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,0,1>();
  }
  else if( (foo[0] == 2) && (foo[1] == 0) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,0,2>();
  }
  else if( (foo[0] == 0) && (foo[1] == 0) && (foo[2] == 0) )
  {
    return radmat::gen_mom<0,0,0>();
  }
  else
  {
    std::cout << __func__ << ": unknown momentum, can " 
      << string_mom(foo) << "  p " << string_mom(p) << std::endl;
  }

}


/////////////////////////////////////////////
  mom_t 
rotate_int_mom(const RotationMatrix_t *R,
    const mom_t &l)
{
  mom_t chk = radmat::gen_mom<0,0,0>();

  for(int i = 0; i < 3; ++i)
  {
    double res(0.); 
    for(int j = 0; j < 3; ++j)
      if ( fabs( (*R)[i+1][j+1] ) > 1e-6 )
        res += (*R)[i+1][j+1] * double(l[j]); 

    chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
  }
  return chk;
} 


/////////////////////////////////////////////
  mom_t 
rotate_z_mom(const RotationMatrix_t *R,
    const double mod_p)
{
  mom_t chk = radmat::gen_mom<0,0,0>();

  for(int i = 0; i < 3; ++i)
  {
    double res(0.); 

    if ( fabs( (*R)[i+1][3] ) > 1e-6 )
      res = (*R)[i+1][3] * mod_p; 

    chk[i] = (res > 0.) ? int( fabs(res) + 0.5 ) : - int( fabs(res) + 0.5 ); 
  }
  return chk;
} 




//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void check_rotations_mom(int argc, char *argv[])
{
  if( argc != 2 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << std::endl; 
  }

  mom_t chk = radmat::gen_mom<0,0,0>();
  int nb(0), ng(0);

  for(int x = -2; x <=2; ++x)
    for(int y = -2; y <=2; ++y)
      for(int z = -2; z <=2; ++z)
      {
        int psq = x*x + y*y + z*z;

        if( psq > 4)
          continue; 

        chk[0] = x;
        chk[1] = y; 
        chk[2] = z; 

        std::string chk_str = string_mom(chk); 
        std::cout << "testing " << chk_str << std::endl;

        RotationMatrix_t *R = radmat::CanonicalRotationEnv::call_factory(chk); 
        mom_t pr = rotate_z_mom(R,sqrt(psq)); 

        if( chk_str == string_mom(pr)  )
        {
          ++ng; 
        }
        else
        {
          std::cout << __func__ << ": Rotation error for " 
            << chk_str << " got " << string_mom(pr) << std::endl;
          ++nb;
          ++ng; 
        }
      } 

  std::cout << ng - nb << " of " 
    << ng << " tests passed " << std::endl;

}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void check_lattice_rotations_mom(int argc, char *argv[])
{
  if( argc != 2 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << std::endl; 
  }

  mom_t chk = radmat::gen_mom<0,0,0>();
  int nb(0), ng(0);

  for(int x = -2; x <=2; ++x)
    for(int y = -2; y <=2; ++y)
      for(int z = -2; z <=2; ++z)
      {
        int psq = x*x + y*y + z*z;

        if( psq > 4)
          continue; 

        chk[0] = x;
        chk[1] = y; 
        chk[2] = z; 

        std::string chk_str = string_mom(chk); 
        std::cout << "testing " << chk_str << std::endl;

        mom_t ref = ref_mom(FF::canonicalOrder(chk)); 
        std::string LG = Hadron::generateLittleGroup(ref); 
        RotationMatrix_t *Rref = radmat::CanonicalRotationEnv::call_factory(LG);  
        RotationMatrix_t *R = radmat::CanonicalLatticeRotationEnv::call_factory(chk); 
        mom_t pr = rotate_z_mom(Rref,sqrt(psq)); 

        if( string_mom(ref) == string_mom(pr) )
        {
          ++ng; 
        }
        else
        {
          std::cout << __func__ << ": Rotation error for " 
            << string_mom(ref) << ": got " 
           << string_mom(pr) << "\nR:\n" << *Rref << std::endl;
          ++nb;
          ++ng; 
        }

        mom_t p = rotate_int_mom(R,pr); 

        if( string_mom(p) == chk_str)
        {
          ++ng;
        }
        else
        {
          std::cout << __func__ << ": rotation error, tried to rotate " 
            << string_mom(pr) << " to " << chk_str << " got " 
            << string_mom(p) << "\nR:\n" << *R << std::endl;
          ++nb;
          ++ng; 
        }

        delete R;
        delete Rref; 
      } 

  std::cout << ng - nb << " of " 
    << ng << " tests passed " << std::endl;

}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void read_momentum(const int st, mom_t &p, char *argv[])
{
  for(int i = st; i < st+3; ++i)
  {
    std::istringstream val(argv[i]);
    val >> p[i-st];
  }
}



void get_lattice_rotation(int argc, char *argv[])
{
  if( argc != 8 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> " << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);

  std::cout << "read m1 " << string_mom(m1) 
    << " m2 " << string_mom(m2) << std::endl;
  std::cout << " returning R_lat * m1 = m2 " << std::endl;

  RotationMatrix_t *r = radmat::generate_frame_transformation(m2,m1);

  if( !!! radmat::check_frame_transformation(r, m1 , m2 ) )
  {
    std::cout << __func__ 
      << ": WARNING -- lattice transformation failed " << std::endl;
  }

  std::cout << "R:" << *r << std::endl;

  std::cout << "det(R) = " << determinant(r) << std::endl;

  delete r; 
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////


void get_frame_rotation(int argc, char *argv[])
{
  if( argc != 14 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <mom_prime1> <mom_prime2>" << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 
  mom_t ma = radmat::gen_mom<0,0,0>(); 
  mom_t mb = radmat::gen_mom<0,0,0>(); 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);
  read_momentum(8,ma,argv);
  read_momentum(11,mb,argv);

  std::cout << "read m1 " << string_mom(m1) << " m2 " << string_mom(m2) 
    << " m_prime1 " << string_mom(ma) << " m_prime2 " << string_mom(mb) << std::endl;

  RotationMatrix_t *rl = radmat::LatticeRotationEnv::get_left_rotation(m1,m2,ma,mb);
  RotationMatrix_t *rr = radmat::LatticeRotationEnv::get_right_rotation(m1,m2,ma,mb);
  std::cout << "Rleft:" << *rl << std::endl;
  std::cout << "Rright:" << *rr << std::endl;

  delete rl; 
  delete rr; 
}


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////



void get_wigner_phase(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <J> " << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 
  int J; 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);
  std::istringstream val(argv[8]);
  val >> J; 


  radmat::LatticeRotationEnv::FrameOrientation_t frame; 

  std::cout << "read m1 " << string_mom(m1) 
    << " m2 " << string_mom(m2) << " J " << J <<std::endl;

  radmat::WignerMatrix_t D; 
  radmat::primitiveEmbededWignerDMatrix *W;  

  if( J == 1 )
    W = new radmat::embededWignerDMatrix<1>();  
  else if (J == 2)
    W = new radmat::embededWignerDMatrix<2>();  
  else if (J == 3)
    W = new radmat::embededWignerDMatrix<3>();  
  else if (J == 4)
    W = new radmat::embededWignerDMatrix<4>();  
  else
  {
    std::cout << "J " << J << "is not supported" << std::endl;
    exit(1); 
  }

  frame = W->get_frame(m1,m2);
  std::cout << "canonical frame " << string_mom(frame.cl) << " " 
    << string_mom(frame.cr) << std::endl;

  D = W->composite_wigner_matrix(m1,m2); 

  std::cout << "\nD:" << D << std::endl;

  radmat::WignerMatrix_t I,Dstar = D;
  W->conjugate(&Dstar); 
  D.lower_index(1); 
  I = radmat::contract(D,Dstar,1,1); 
  W->clean_up(I); 

  std::cout << "D Ddag = " << I << std::endl;


  delete W; 
}




//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_left_triad_wigner(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <J> " << std::endl; 
    exit(1);
  }
  
  mom_t l = radmat::gen_mom<0,0,0>(); 
  mom_t r = radmat::gen_mom<0,0,0>(); 
  int J; 

  read_momentum(2,l,argv);
  read_momentum(5,r,argv);
  std::istringstream val(argv[8]);
  val >> J; 

  radmat::DMatrixManager Wig; 
  radmat::RotationMatrix_t * Rtriad; 
  Rtriad = Wig.triad_rotation_matrix(l,r);
  radmat::WignerMatrix_t *D,*Dtriad; 
  D = Wig.left_wigner_matrix(Rtriad,l,r,J,true); 

  clean_up_rot_mat(Rtriad);
  Wig.clean(D); 

  std::cout << "l " << string_mom(l) << " r " << string_mom(r)
    << "\nD_left: " << *D << "\nRtriad:" << *Rtriad << std::endl;  

  delete D; 
  delete Rtriad; 
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_right_triad_wigner(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom1> <mom2> <J> " << std::endl; 
    exit(1);
  }
  
  mom_t l = radmat::gen_mom<0,0,0>(); 
  mom_t r = radmat::gen_mom<0,0,0>(); 
  int J; 

  read_momentum(2,l,argv);
  read_momentum(5,r,argv);
  std::istringstream val(argv[8]);
  val >> J; 


  radmat::DMatrixManager Wig; 
  radmat::RotationMatrix_t * Rtriad; 
  Rtriad = Wig.triad_rotation_matrix(l,r);
  radmat::WignerMatrix_t *D,*Dtriad; 
  D = Wig.right_wigner_matrix(Rtriad,l,r,J,true); 

  clean_up_rot_mat(Rtriad);
  Wig.clean(D); 

  std::cout << "l " << string_mom(l) << " r " << string_mom(r)
    << "\nD_left: " << *D << "\nRtriad:" << *Rtriad << std::endl;  

  delete D; 
  delete Rtriad; 
}





//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_wigner_matrix(int argc, char *argv[])
{
  if( argc != 6 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom> <J> " << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  int J; 

  std::istringstream val(argv[5]);
  val >> J; 

  read_momentum(2,m1,argv);

  std::cout << "read p " << string_mom(m1) 
    << " J " << J << std::endl;
  std::cout << " Wigner_D(p,J) " << std::endl;

  radmat::WignerMatrix_t *W = radmat::WignerDMatrixEnv::call_factory(m1,J);

  std::cout << "W:" << *W << std::endl;

  delete W; 
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

void get_prod_wigner_matrix(int argc, char *argv[])
{
  if( argc != 9 )
  {
    std::cout << "error usage: test_rotations [" << __func__ << "] "
      << "<mom> <mom> <J> " << std::endl; 
    exit(1);
  }
  
  mom_t m1 = radmat::gen_mom<0,0,0>(); 
  mom_t m2 = radmat::gen_mom<0,0,0>(); 
  int J; 

  std::istringstream val(argv[8]);
  val >> J; 

  read_momentum(2,m1,argv);
  read_momentum(5,m2,argv);

  std::cout << "read p " << string_mom(m1) 
    << " J " << J << std::endl;
  std::cout << " Wigner_D(p,J) * Wigner_D(p2,J) " << std::endl;

  radmat::WignerMatrix_t *W = radmat::WignerDMatrixEnv::call_factory(m1,J);
  radmat::WignerMatrix_t *W2 = radmat::WignerDMatrixEnv::call_factory(m2,J);

  W->lower_index(1); 
  radmat::WignerMatrix_t WW = radmat::contract(*W,*W2,1,0); 

  std::cout << "W:" << WW << std::endl;

  delete W; 
  delete W2; 
}



//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

typedef void (*fptr)(int argc , char *argv[]); 
std::map<std::string , fptr > options;

void insert_op(const std::string &s, const fptr &f)
{
  options.insert(std::pair<std::string,fptr>(s,f)); 
}

void init_options(void)
{
  insert_op("get_left_triad_wigner",&get_left_triad_wigner);
  insert_op("get_right_triad_wigner",&get_right_triad_wigner);
  insert_op("check_rotations_mom",&check_rotations_mom); 
  insert_op("check_lattice_rotations_mom",&check_lattice_rotations_mom); 
  insert_op("get_lattice_rotation",&get_lattice_rotation); 
  insert_op("get_frame_rotation",&get_frame_rotation); 
  insert_op("get_wigner_matrix",&get_wigner_matrix); 
  insert_op("get_prod_wigner_matrix",&get_prod_wigner_matrix); 
  insert_op("get_wigner_phase",&get_wigner_phase); 
}


// pick appropriate function and pass on command line inputs 
void do_work(std::string &op, int argc,char *argv[])
{
  init_options(); 

  if(options.find(op) == options.end())
  {
    std::cerr << " unrecognized op " << op 
      << " more intelligent choices are " << std::endl; 
    std::map<std::string , fptr>::const_iterator it; 
    for(it = options.begin(); it != options.end(); ++it)
      std::cerr << it->first << std::endl; 
    exit(1); 
  }

  fptr foo = options[op];

  foo(argc,argv); 
}


// main program wrapper
int main(int argc, char *argv[])
{
  radmat::AllFactoryEnv::registerAll(); 

  // we will always have at least 2 , radmat_util operation_with_no_inputs
  if(argc < 2)
  {
    std::cerr << "usage: test_rotations : <operation> <op inputs ...> " << std::endl;
    exit(1); 
  }

  std::string op;
  std::istringstream opi(argv[1]); 
  opi >> op; 

  do_work(op,argc,argv); 

  return 0;
}










