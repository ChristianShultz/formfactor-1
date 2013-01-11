/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : mass_overlap_db_create.cc

 * Purpose :

 * Creation Date : 08-01-2013

 * Last Modified : Tue Jan  8 13:08:13 2013

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
#include "AllConfStoreDB.h"
#include "ensem/ensem.h"
#include "radmat/load_data/radmat_overlap_key_val_db.h"
#include "hadron/hadron_npart_irrep.h"
#include "semble/semble_meta.h"


using namespace radmat;
using namespace Hadron;
using namespace ENSEM;
using namespace ADATXML;


struct XML_input_t
{
  std::string mass_file;
  std::string overlap_file;
  int ncfg;
  bool resize;
  bool isProjected;
  std::string dbname;
  std::string pid; 
  KeyHadronNPartIrrep_t redstar; 
};


namespace
{

  template<typename T>
    void doXMLRead(ADATXML::XMLReader &ptop, const std::string &path, T &place, const char * f)
    {
      if(ptop.count(path) > 0)
        read(ptop,path,place);
      else
      {
        std::cerr << __PRETTY_FUNCTION__ << ": Error, called by " << f << " trying to read path, " << path
          << ", path was empty, exiting" << std::endl;
        exit(1);
      }
    }
} // namespace anonomyous 


void read(ADATXML::XMLReader &xml, const std::string &path, XML_input_t &prop)
{
  ADATXML::XMLReader ptop(xml,path);
  doXMLRead(ptop,"mass_file",prop.mass_file,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"overlap_file",prop.overlap_file,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"ncfg",prop.ncfg,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"resize",prop.resize,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"isProjected",prop.isProjected,__PRETTY_FUNCTION__); 
  doXMLRead(ptop,"dbname",prop.dbname,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"pid",prop.pid,__PRETTY_FUNCTION__);
  doXMLRead(ptop,"redstar",prop.redstar,__PRETTY_FUNCTION__);
}



struct dbInterface
{
  typedef RadmatExtendedKeyHadronNPartIrrep_t K;
  typedef RadmatMassOverlapData_t D;
  typedef ADATIO::SerialDBKey<K> SK;
  typedef ADATIO::SerialDBData<D> SD;

  dbInterface(void); // hide ctor
  dbInterface(const XML_input_t &xml) : m_xml(xml) {} 
  ~dbInterface(void) {if (m_db); delete m_db; m_db = NULL;}

  bool alloc(void)
  {
    try
    {
      m_db = new FILEDB::AllConfStoreDB<SK,SD>(std::vector<int>(1,1));
    }
    catch(...)
    {
      std::cerr << __func__ << ": ERROR: couldn't alloc a database" << std::endl;
      delete m_db;
      m_db = NULL;
      exit(1);
    }

    if(m_db->open(m_xml.dbname,  O_RDWR | O_TRUNC | O_CREAT, 0664) != 0)
    {
      std::cerr << __func__ << ": error opening dbase= " << m_xml.dbname << std::endl;
      exit(1);
    }
    return true;
  }


  void insert(void)
  {
    if(!!!m_db)
      if(!!!alloc())
      {
        std::cerr << __func__ << ": ERROR: something bad happened" << std::endl;
        exit(1);
      }

    SK key;
    SD data;

    key.key() = K(m_xml.pid,m_xml.redstar);
    std::stringstream E,Z;
    ENSEM::read(m_xml.mass_file,data.data().E());
    ENSEM::read(m_xml.overlap_file,data.data().Z());


    if(data.data().E().size() < m_xml.ncfg)
    {
      std::cout << "Error: don't have enough configurations on the masses/overlaps to do " 
        << "the problem.. exiting" << std::endl;
      exit(1);
    }

    // allow for some sort of hackey resampling if we didn't generate enought 3pt statistics
    if(data.data().E().size() > m_xml.ncfg)
      if(m_xml.resize)
      {
        ENSEM::EnsemReal tmpE,tmpZ,tmpEE,tmpZZ;
        
        tmpE = ENSEM::rescaleEnsemDown(data.data().E());
        tmpZ = ENSEM::rescaleEnsemDown(data.data().Z());
        tmpEE.resize(m_xml.ncfg);
        tmpZZ.resize(m_xml.ncfg);

        for(int i = 0; i < m_xml.ncfg; ++i)
        {
          tmpEE.elem(i) = tmpE.elem(i);
          tmpZZ.elem(i) = tmpZ.elem(i);
        }

        data.data().E() = ENSEM::rescaleEnsemUp(tmpEE); 
        data.data().Z() = ENSEM::rescaleEnsemUp(tmpZZ); 
      }

    if(m_xml.isProjected)
      data.data().Z() = SEMBLE::toScalar(1.); 

    if(m_db->insert(key,std::vector<SD>(1,data)) != 0)
    {
      std::cerr << __func__ << ": could not insert key \n" << key.key() << std::endl;
      exit(1);
    }
  }

  XML_input_t m_xml; 
  FILEDB::AllConfStoreDB< SK , SD > *m_db;
};


int main(int argc , char *argv[])
{
  if(argc != 2)
  {
    std::cerr << "usage: " << argv[0] << ": <input.xml>" << std::endl;
    exit(1);
  }

  XML_input_t ini;
  std::string xmlinifile;
  std::stringstream val(argv[1]);
  val >> xmlinifile;

  try
  {
    ADATXML::XMLReader xml(xmlinifile);
    read(xml,"/Params",ini);
  }
  catch(const std::string &e)
  {
    std::cerr << __func__ << ": ERROR: can't read xmlinifile (" << xmlinifile << "): " << e << std::endl;
    exit(1);
  }


  dbInterface foo(ini); 
  
  if(!!!foo.alloc())
  {
    std::cout << "something bad happened" << std::endl;
    exit(1);
  }

  foo.insert();


  return 0;
}

