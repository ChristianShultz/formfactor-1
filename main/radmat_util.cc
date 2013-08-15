/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : radmat.cc

 * Purpose :

 * Creation Date : 25-02-2013

 * Last Modified : Tue Jun 25 17:04:07 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include "radmat/driver/radmat_driver.h"




void gen_xml(int argc , char *argv[] )
{
  if(argc != 4)
  {
    std::cerr << "error: usage: radmat_util: gen_xml <xmlinifile> <mode> " << std::endl;
    exit(1); 
  }

  std::istringstream val(argv[2]); 
  std::string ini; 
  val >> ini; 

  std::istringstream val2(argv[3]); 
  std::string mode; 
  val2 >> mode; 

  radmat::RadmatDriver d; 
  d.xml_handler(ini,mode); 
}


void nuke_graphs(int argc , char *argv[] )
{

  if(argc != 5)
  {
    std::cerr << "error: usage: radmat_util: nuke_graphs" 
      <<" <smeared_xmlinifile> <smeared_graph.sdb> <output.xml>" << std::endl; 

    std::cerr << " this takes an xml ini with a smeared insertion and the graph db \n" 
      << " that redstar_gen_graph spits out and comes up with an xml list of the \n" 
      << " disconnected insertions that need to be nuked" << std::endl;
    exit(1); 
  }


  std::istringstream val1(argv[2]);
  std::istringstream val2(argv[3]); 
  std::istringstream val3(argv[4]); 
  std::string ini,graph,xmlout;
  val1 >> ini; 
  val2 >> graph; 
  val3 >> xmlout; 

  radmat::RadmatDriver d; 
  d.nuke_graph(ini,graph,xmlout); 

}

void stub_xml(int argc , char *argv[] )
{
  if(argc != 3)
  {
    std::cerr << "error: usage: radmat_util: stub_xml <xmlinifile> " << std::endl;
    exit(1); 
  }

  std::istringstream val(argv[2]); 
  std::string ini; 
  val >> ini; 

  radmat::RadmatDriver d; 
  d.build_stub_xml(ini); 

}



// typedef 
typedef void (*fptr)(int argc , char *argv[]) ; 

// a map of operation names and function pointers
std::map<std::string , fptr> options; 

// init the map 
void init_options(void)
{
  options.insert(std::pair<std::string,fptr>("gen_xml",&gen_xml)); 
  options.insert(std::pair<std::string,fptr>("nuke_graphs",&nuke_graphs));
  options.insert(std::pair<std::string,fptr>("stub_xml",&stub_xml));  
}

// pick appropriate function and pass on command line inputs 
void do_work(std::string &op, int argc,char *argv[])
{
  init_options(); 

  if(options.find(op) == options.end())
  {
    std::cerr << " unrecognized op " << op << " intelligent choices were " << std::endl; 
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
  // we will always have at least 2 , radmat_util operation_with_no_inputs
  if(argc < 2)
  {
    std::cerr << "usage: radmat_util : <operation> <op inputs ...> " << std::endl;
    exit(1); 
  }

  std::string op;
  std::istringstream opi(argv[1]); 
  opi >> op; 

  do_work(op,argc,argv); 


  return 0;
}
