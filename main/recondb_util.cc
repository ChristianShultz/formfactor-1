/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : recondb_util.cc

 * Purpose :

 * Creation Date : 05-10-2012

 * Last Modified : Tue Dec 11 09:39:14 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/






#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>
#include <vector>
#include <utility>
#include <map>
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"
#include "semble_fit_ini_xml.h"
#include "hadron/hadron_npart_irrep.h"
#include "ensem/ensem.h"
#include "semble/semble_key_val_db.h"
#include "semble/semble_meta.h"
#include "AllConfStoreDB.h"


typedef SEMBLE::SembleExtendedKeyHadronNPartIrrep_t K;
typedef SEMBLE::SembleMassOverlapData_t D;
typedef ADATIO::SerialDBKey<K> SK;
typedef ADATIO::SerialDBData<D> SD;


template<typename KEY, typename DATA> 
  typename std::vector< ADATIO::SerialDBKey<KEY> >
keys(typename FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > &db)
{
  typename std::vector< ADATIO::SerialDBKey<KEY> > skeys;
  db.keys(skeys);
  return skeys;
}

template<typename KEY, typename DATA> 
  void
print_keys(typename FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > &db)
{
  typename std::vector< ADATIO::SerialDBKey<KEY> > m_keys = keys(db);
  typename std::vector< ADATIO::SerialDBKey<KEY> >::const_iterator it;

  for(it = m_keys.begin(); it != m_keys.end(); ++it)
    std::cout << it->key() << std::endl;
}

template<typename KEY, typename DATA>
  void
print_keys_xml(typename FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > &db, const std::string &xmlfile)
{
  typename std::vector< ADATIO::SerialDBKey<KEY> > m_keys = keys(db);
  ADATXML::Array<KEY> xml_keys; 
  xml_keys.resize(m_keys.size());

  for(int i = 0; i < m_keys.size(); ++i)
    xml_keys[i] = m_keys[i].key(); 

  ADATXML::XMLBufferWriter buff; 
  write(buff,"Keys",xml_keys);
  std::ofstream out(xmlfile.c_str());
  buff.print(out);
  out.close();
}


void blow_up(const char *c)
{
  std::cout << __func__ << ": error: called by " << c << std::endl;
  exit(1);
}

  template<typename KEY, typename DATA>
int work_handler_2par(const std::string &dbname, const std::string &op)
{

  FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > db;
  if(db.open(dbname, O_RDONLY , 0400) != 0)
  {
    std::cerr << __func__ << ": error opening database " << dbname << std::endl;
    exit(1);
  }

  if(op == "keys")
    print_keys(db);
  else
    blow_up(__PRETTY_FUNCTION__); 

  return 0;
}

  template<typename KEY, typename DATA>
int work_handler_3par(const std::string &dbname, const std::string &op, const std::string &par)
{

  FILEDB::AllConfStoreDB< ADATIO::SerialDBKey<KEY> , ADATIO::SerialDBData<DATA> > db;
  if(db.open(dbname, O_RDONLY , 0400) != 0)
  {
    std::cerr << __func__ << ": error opening database " << dbname << std::endl;
    exit(1);
  }

  if(op == "keys_xml")
    print_keys_xml(db,par);
  else
    blow_up(__PRETTY_FUNCTION__); 

  return 0;
}



int main(int argc , char *argv[])
{

  if(argc <= 3)
  {
    std::cerr << "usage; " << argv[0] << ": <database> <operation> [..<operation options>]" << std::endl;
    exit(1);
  }

  // Read the dbname
  std::string dbname;
  std::istringstream dbnamestream(argv[1]);
  dbnamestream >> dbname;

  std::string op;
  std::istringstream opstream(argv[2]);
  opstream >> op;

  if(argc == 3)
    return work_handler_2par<K,D>(dbname,op);

  std::string par1;
  std::istringstream par1stream(argv[3]);
  par1stream >> par1;

  if(argc == 4)
  {
    return work_handler_3par<K,D>(dbname,op,par1); 
  }


  blow_up(__PRETTY_FUNCTION__); 
} 





