/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : read_redstar.cc

 * Purpose :

 * Creation Date : 14-08-2012

 * Last Modified : Mon Aug 20 15:07:56 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/




#include "io/adat_xmlio.h"
#include "io/key_val_db.h"

#include "DBString.h"

#include "hadron/prop_elem_type.h"
#include "hadron/genprop_elem_type.h"
#include "hadron/meson_elem_type.h"
#include "hadron/glueball_elem_type.h"
#include "hadron/baryon_elem_type.h"
#include "formfac/hadron_1pt_corr.h"
#include "formfac/hadron_2pt_corr.h"
#include "formfac/hadron_3pt_corr.h"

#include "hadron/hadron_node.h"
#include "hadron/hadron_npart_npt_conn_graph.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "hadron/hadron_timeslice.h"

#include "radmat/utils/splash.h"
#include "radmat/load_data/load_database.h"


#include <string>
#include <sstream>
#include <exception>


using namespace std;
using namespace ENSEM;
using namespace ADATXML;
using namespace ADATIO;
using namespace FILEDB;
using namespace FF;
using namespace Hadron;


typedef KeyHadronNPartNPtCorr_t K;
typedef EnsemVectorComplex V;

void simpleWorld(const std::string &dbase)
{

  ConfDataStoreDB< SerialDBKey<K>,  SerialDBData<V> > database;
  database.open(dbase, O_RDONLY, 0400);

  std::vector< SerialDBKey<K> > keys;
  database.keys(keys);

  for(unsigned int i=0; i < keys.size(); ++i)
  {
    std::cout << keys[i].key();
  }
}



std::vector<SerialDBKey<K> > getSerialKeys(const std::string &dbase)
{
  ConfDataStoreDB< SerialDBKey<K>,  SerialDBData<V> > database;
  database.open(dbase, O_RDONLY, 0400);

  std::vector< SerialDBKey<K> > keys;
  database.keys(keys);

  return keys;
}

std::vector<K> getKeys(const std::string &dbase)
{
  std::vector<SerialDBKey<K> > serial_keys = getSerialKeys(dbase);
  std::vector<SerialDBKey<K> >::const_iterator it;
  std::vector<K> keys;

  for(it = serial_keys.begin(); it != serial_keys.end(); ++it)
    keys.push_back(it->key());


  return keys;
}

void lessSimpleWorld(const std::string &dbase)
{

  radmat::LoadFromDatabase<K,V> m_database_loader(dbase);
  std::vector<K> keys = getKeys(dbase); 
  K key = keys.at(0);
  V val = m_database_loader.getFromKey(key);
  std::vector<V> vals = m_database_loader.getFromKey(keys);

  for(unsigned int i = 0; i < vals.size(); ++i)
  {
    std::cout << keys[i] << std::endl;
    std::cout << ENSEM::mean(vals[i]) << std::endl;
  }


}


int main(int argc , char* argv[])
{

  if(argc != 2)
  {
    SPLASH("usage: read_redstar : <dbfile.sdb>");
    exit(1);
  }

  std::istringstream fname(argv[1]);
  std::string name;
  fname >> name;

  std::cout << "Loading database: " << name << std::endl;

  try
  {
    simpleWorld(name);
  }
  catch(...)
  {
    SPLASH("An exception occurred in this context");
    exit(1);
  }


  try
  {
    lessSimpleWorld(name);
  }
  catch(...)
  {
    SPLASH("An exception occurred in this context");
    exit(1);
  }


  return 0;
}
