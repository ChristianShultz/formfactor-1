/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : print_twoquark_ops.cc

 * Purpose :

 * Creation Date : 23-11-2013

 * Last Modified : Sat 23 Nov 2013 01:58:15 PM EST

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "hadron/twoquark_dirac_ops.h"
#include "hadron/twoquark_displace.h"
#include "ensem/ensem.h"
#include "semble/semble_meta.h"
#include "itpp/itbase.h"
#include <map>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>
#include <utility>

// the junk in adat/lib/hadron
typedef Hadron::TheTwoQuarkDiracSpinOpsRegInfoFactory RegMap; 
typedef Hadron::TheTwoQuarkDiracSpinOpsFactory TwoQuarkFoundry; 


// a hacky thingy
struct KeyRetriever 
{
  template<typename T>
    typename T::first_type operator() (T KeyValPair) const
    {
      return KeyValPair.first; 
    }
};

// what are we trying to do 
  std::pair<std::string,int>
usage(int argc, char *argv[] ) 
{
  if ( argc != 3 ) 
  {
    std::cerr << "error usage: print_twoquark_ops <pattern> <row> " << std::endl; 
    exit(1337); 
  }

  std::istringstream val(argv[1]), val2(argv[2]);
  std::string pattern; 
  int row; 
  val >> pattern; 
  val2 >> row; 
  std::cout << "Matching on " << pattern 
    << " row " << row << std::endl ; 


  return std::pair<std::string,int>(pattern,row);
}

  itpp::Mat<std::complex<double> > 
convertSpinMatrix(const ENSEM::SpinMatrix &m)
{
  itpp::Mat<std::complex<double> > spin_mat(4,4);
  spin_mat.zeros(); 

  for(int i =0; i < 4; ++i)
    for(int j = 0; j<4; ++j) 
      spin_mat(i,j) = SEMBLE::toScalar(ENSEM::peekSpin(m,i,j)); 

  return itpp::round_to_zero(spin_mat,1e-6); 
}

std::string do_print(const std::complex<double> &c)
{
  double r = c.real(); 
  double i = c.imag(); 

  std::stringstream ss; 

  if ( r != 0.) 
  {
    ss << r << " ";
    if(i != 0) 
    {
      if ( i > 0) 
        ss << "+ " << i << "i ";
      else 
        ss << "- " << fabs(i) << "i " ; 
    }
  }
  else if (i != 0) 
  {
    if ( fabs(i) == 1. )
    {
      if ( i > 0) 
        ss << "i ";
      else 
        ss << "- " << "i " ; 
    } 
    else
      ss << i << "i "; 
  }
  else
  {
    ss << 0 << " "  ; 
  }

  return ss.str(); 
}


void printSpinMatrix(const ENSEM::SpinMatrix &m)
{
  itpp::Mat<std::complex<double> > spin_matrix = convertSpinMatrix(m); 
  std::cout << "\\begin{bmatrix}" << std::endl; 

  for(int i =0; i < 4; ++i)
  {
    for(int j = 0; j<3; ++j) 
      std::cout << do_print(spin_matrix(i,j)) << " & "; 

    std::cout << do_print(spin_matrix(i,3))  << "\\\\" << std::endl;  
  }
  std::cout << "\\end{bmatrix}\n" << std::endl; 
}


// do work 
int main(int argc, char *argv[])
{

  // get pattern from cmd line  
  std::pair<std::string,int> u = usage(argc,argv); 
  std::string pattern = u.first; 
  int row = u.second; 

  // register ops
  Hadron::TwoQuarkDiracSpinOpsEnv::registerAll(); 

  // keys and iterator
  std::vector<std::string> reg_keys; 
  std::vector<std::string>::const_iterator  reg_keys_it; 

  // pull the keys out 
  std::transform(RegMap::Instance().begin(),
      RegMap::Instance().end(),
      std::back_inserter(reg_keys),
      KeyRetriever()); 

  // print them
  for(reg_keys_it = reg_keys.begin(); reg_keys_it != reg_keys.end(); ++reg_keys_it)
    if( reg_keys_it->find(pattern) != std::string::npos)
    {
      std::cout << "matched -> " <<  *reg_keys_it << std::endl; 

      Hadron::TwoQuarkDiracSpinOpsRegInfo_t reg = RegMap::Instance().find(*reg_keys_it)->second;
      Hadron::TwoQuarkDiracSpinOps* op = TwoQuarkFoundry::Instance().createObject(reg.id); 

      Hadron::MapTwoQuarkExpr_t expr =  (*op)(Hadron::initOne(),row); 

      std::vector<Hadron::KeyTwoQuarkExpr_t> expr_keys = expr.keys(); 
      std::vector<Hadron::KeyTwoQuarkExpr_t>::const_iterator k;

      for( k = expr_keys.begin(); k != expr_keys.end(); ++k)
        printSpinMatrix(expr[*k]);
    }


  // now make the dictionary 
  for(int i = 0; i < 4; ++i) 
  {
    Hadron::MapTwoQuarkExpr_t expr = ENSEM::GammaDP(1 << i) * Hadron::initOne(); 
    std::vector<Hadron::KeyTwoQuarkExpr_t> expr_keys = expr.keys(); 
    std::vector<Hadron::KeyTwoQuarkExpr_t>::const_iterator k;
    std::cout << "GammaDP(" << i << "):" << std::endl; 
    for( k = expr_keys.begin(); k != expr_keys.end(); ++k)
      printSpinMatrix(expr[*k]);
  }

  return 0; 
}

