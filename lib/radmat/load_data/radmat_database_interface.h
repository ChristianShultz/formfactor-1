#ifndef RADMAT_DATABASE_INTERFACE_H
#define RADMAT_DATABASE_INTERFACE_H


#include <vector>
#include <exception>
#include "ensem/ensem.h"
#include "adat/handle.h"
#include "AllConfStoreDB.h"
#include "io/adat_io.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"


namespace radmat
{

  struct dbProp_t
  {
    std::string dbname;
    std::string badlist;
  };


  std::string toString(const dbProp_t &);
  std::ostream& operator<<(std::ostream&,const dbProp_t&);
  void read(ADATXML::XMLReader &xml, const std::string &path, dbProp_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const dbProp_t &);


  struct radmatDBProp_t
  {
    dbProp_t threePointDatabase;
    dbProp_t normalizationDatabase;
  };

  std::string toString(const radmatDBProp_t &);
  std::ostream& operator<<(std::ostream&, const radmatDBProp_t &);
  void read(ADATXML::XMLReader &xml, const std::string &path, radmatDBProp_t &);
  void write(ADATXML::XMLWriter &xml, const std::string &path, const radmatDBProp_t &);


  // some templated structure to deal with all of the database nastyness 
  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    struct radmatAllConfDatabaseInterface
    {
      typedef typename ADATIO::SerialDBKey<CORRKEY> S_C_KEY;
      typedef typename ENSEM::EnsemScalar<CORRDATA>::Type_t ESCALARDATA;
      typedef typename ADATIO::SerialDBData<ESCALARDATA> S_C_DATA;
      typedef typename ADATIO::SerialDBKey<NORMKEY> S_N_KEY;
      typedef typename ADATIO::SerialDBData<NORMDATA> S_N_DATA; 

      radmatAllConfDatabaseInterface(void); // hide ctor
      radmatAllConfDatabaseInterface(const radmatDBProp_t &);


      void alloc(void);
      bool exists(const CORRKEY &) const;
      bool exists(const NORMKEY &) const;
      CORRDATA fetch(const CORRKEY &) const;
      NORMDATA fetch(const NORMKEY &) const; 


      ADAT::Handle<FILEDB::AllConfStoreDB<S_C_KEY,S_C_DATA> > m_corr_db;
      ADAT::Handle<FILEDB::AllConfStoreDB<S_N_KEY,S_N_DATA> > m_norm_db;  
      radmatDBProp_t db_props; 
    };

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA >
    radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::radmatAllConfDatabaseInterface(const radmatDBProp_t &prop)
    : db_props(prop)
    {
      alloc();
    }


  // set up some handles and check that we can open it
  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    void radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::alloc(void)
    {

      m_corr_db = ADAT::Handle<FILEDB::AllConfStoreDB<S_C_KEY,S_C_DATA> > (new FILEDB::AllConfStoreDB<S_C_KEY,S_C_DATA>() );
      m_norm_db = ADAT::Handle<FILEDB::AllConfStoreDB<S_N_KEY,S_N_DATA> > (new FILEDB::AllConfStoreDB<S_N_KEY,S_N_DATA>() );

      if(m_corr_db->open(db_props.threePointDatabase.dbname, O_RDWR | O_TRUNC | O_CREAT, 0664) != 0)
      {
        std::cerr << __func__ << ": error opening dbase= " << db_props.threePointDatabase.dbname << std::endl;
        exit(1);
      }


      if(m_norm_db->open(db_props.normalizationDatabase.dbname, O_RDWR | O_TRUNC | O_CREAT, 0664) != 0)
      {
        std::cerr << __func__ << ": error opening dbase= " << db_props.normalizationDatabase.dbname << std::endl;
        exit(1);
      }
    }

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    bool radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::exists(const CORRKEY &k) const
    {
      S_C_KEY key;
      key.key() = k;
      return m_corr_db->exist(key);     
    }

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    bool radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::exists(const NORMKEY &k) const
    {
      S_N_KEY key;
      key.key() = k;
      return m_norm_db->exist(key);
    }

  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    CORRDATA  radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::fetch(const CORRKEY &k) const
    {
      CORRDATA eval;
      try
      {
        S_C_KEY key;
        key.key() = k;
        std::vector<S_C_DATA> vals;
        int ret(0);
        if ((ret = m_corr_db->get(key, vals)) != 0)
        {
          std::cerr << __func__ << ": key not found\n" << k;
          exit(1);
        }

        eval.resize(vals.size());
        eval.resizeObs(vals[0].data().numElem());

        for(unsigned int i=0; i < vals.size(); ++i)
        {
          ESCALARDATA sval = vals[i].data();
          ENSEM::pokeEnsem(eval, sval, i);
        }

      }
      catch(const std::string& e) 
      {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
      }
      catch(std::exception& e) 
      {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
      }

      return eval;
    }


  template<typename CORRKEY, typename CORRDATA, typename NORMKEY, typename NORMDATA>
    NORMDATA  radmatAllConfDatabaseInterface<CORRKEY,CORRDATA,NORMKEY,NORMDATA>::fetch(const NORMKEY &k) const
    {
      NORMDATA eval;
      try
      {
        S_N_KEY key;
        key.key() = k;
        std::vector<S_N_DATA> val; // hacky thingy we did
        int ret(0);
        if ((ret = m_norm_db->get(key, val)) != 0)
        {
          std::cerr << __func__ << ": key not found\n" << k;
          exit(1);
        }

        eval = val[0].data(); // if this worked right we stored the entire thing in a single ensemble slot.. 

      }
      catch(const std::string& e) 
      {
        std::cerr << __func__ << ": Caught Exception: " << e << std::endl;
        exit(1);
      }
      catch(std::exception& e) 
      {
        std::cerr << __func__ << ": Caught standard library exception: " << e.what() << std::endl;
        exit(1);
      }

      return eval;
    } 





} // namespace radmat













#endif /* RADMAT_DATABASE_INTERFACE_H */
