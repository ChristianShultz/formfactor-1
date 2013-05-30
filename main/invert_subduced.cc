/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

* File Name : invert_subduced.cc

* Purpose :

* Creation Date : 03-05-2013

* Last Modified : Mon May  6 16:32:45 2013

* Created By : shultz

_._._._._._._._._._._._._._._._._._._._._.*/



#include "radmat/load_data/invert_subduction.h"
#include "semble/semble_meta.h"

#include <string>
#include <sstream>





using namespace std;
using namespace radmat; 


std::map<std::string, bool> pars; 

void init_pars(void)
{
  pars.insert(std::pair<std::string,bool>("pos",true)); 
  pars.insert(std::pair<std::string,bool>("neg",false)); 
}

bool read_par(const std::string &par)
{
  std::map<std::string, bool>::const_iterator it; 
  it = pars.find(par); 

  if(it == pars.end())
  {
    std::cerr << "unrecognized parity, use" << std::endl;
    for(it = pars.begin(); it != pars.end(); ++it)
      std::cerr << it->first << std::endl;

    exit(1); 
  }

  return it->second; 
}


ContinuumBosonExprPrimitive usage(int argc, char *argv[])
{

  if(argc != 5)
  {
    std::cerr << ": usage invert <J> <pos,neg> <H> <group (Oh,D2,D3,D4)>" << std::endl;
    exit(1);  
  }

  init_pars(); 

  int J;
  bool P;
  int H;
  std::string par; 
  std::string group; 

  {std::stringstream val(argv[1]); val >> J;}
  {std::stringstream val(argv[2]); val >> par; P = read_par(par);}
  {std::stringstream val(argv[3]); val >> H;}
  {std::stringstream val(argv[4]); val >> group;} 


    return ContinuumBosonExprPrimitive(J,P,H,group); 
}



int main(int argc , char *argv[])
{
  
  ContinuumBosonExprPrimitive fred = usage(argc,argv); 
  
  ListLatticeIrrepExpr_t freddy = invertSubduction(fred); 

  ListLatticeIrrepExpr_t::const_iterator it; 
  for(it = freddy.begin(); it != freddy.end(); ++it)
    std::cout << SEMBLE::toScalar(it->m_coeff) << " X " << it->m_obj << std::endl;

  return 0;
}

