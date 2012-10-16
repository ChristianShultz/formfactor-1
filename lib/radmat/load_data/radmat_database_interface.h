#ifndef RADMAT_DATABASE_INTERFACE_H
#define RADMAT_DATABASE_INTERFACE_H


#include <string>
#include <vector>

#include "AllConfStoreDB.h"
#include "io/key_val_db.h"
#include "ensem/ensem.h"

#include "radmat/utils/pow2assert.h"

#include <iostream>
#include <exception>



namespace radmat
{

  template<typename KEY, typename DATA>
    struct RadmatDBInterface
    {
      //! nested to hold path to db, and a key file we may use
      struct DBParams
      { 
        std::string dbFile;  //! path to the database, member of xml input
        std::string keyFile; //! path to a key file, member of xml input, can be blank
        bool useKeyFile;     //! do we want to use the key file 
      };


      //! Default constructor
      RadmatDBInterface(void)
        : m_haveparams(false) 
      {}

      //! constructor taking the parameters we want to use
      RadmatDBInterface(const DBParams &pars)
        : m_params(pars) , m_haveparams(true)
      {}

      //! load or update parameters
      void loadParams(const DBParams &pars) 
      {
        m_params(pars);
        m_haveparams(true);
      }

      //! get the keys, eiter from the file or all keys in the db
      std::vector<KEY> getKeys(void) const;

      //! get the data corresponding to a key one at a time (slow)
      DATA getEnsem(const KEY &k) const; 

      //! get lots of data at once
      std::vector<DATA> getEnsem(const std::vector<KEY> &k) const;

      DBParams m_params;  //! struct holding params from xml input
      bool m_haveparams;  //! internal check variable

    };


  // get the set of keys we will be working with
  template<typename K, typename D>
    std::vector<K> RadmatDBInterface<K,D>::getKeys(void) const
    {
      POW2_ASSERT(m_haveparams);

            typedef typename ENSEM::EnsemScalar<D>::Type_t SD;


      std::vector<K> keys;

      // TODO
      if(m_params.useKeyFile)
      {
        std::cerr << __func__ << "Christian is a moron and forgot to implement this.." 
          <<" either fix it yourself or bother him" << std::endl;
        exit(1);
      }
      else
      {
        FILEDB::AllConfStoreDB<ADATIO::SerialDBKey<K> , ADATIO::SerialDBData<SD> > database;
        if (database.open(m_params.dbFile, O_RDONLY, 0400) != 0)
        {
          std::cerr << __func__ << ": error opening dbase= " << m_params.dbFile << std::endl;
          exit(1);
        }

        std::vector< ADATIO::SerialDBKey<K> > serial_keys;
        database.keys(serial_keys);

        const unsigned int sz = serial_keys.size();

        keys.resize(sz);    

        for(unsigned int index = 0; index < sz; ++index)
          keys[index] = serial_keys[index].key();
      }

      return keys;
    }


  // get the value of a single key
  template<typename K, typename D>
    D RadmatDBInterface<K,D>::getEnsem(const K &k) const
    {
      typedef typename ENSEM::EnsemScalar<D>::Type_t SD;

      POW2_ASSERT(m_haveparams);

      FILEDB::AllConfStoreDB<ADATIO::SerialDBKey<K> , ADATIO::SerialDBData<SD> > database;
      if (database.open(m_params.dbFile, O_RDONLY, 0400) != 0)
      {
        std::cerr << __func__ << ": error opening dbase= " << m_params.dbFile << std::endl;
        exit(1);
      }

      D eval;

      try
      {
        ADATIO::SerialDBKey<K> key;
        key.key() = k;

        std::vector< ADATIO::SerialDBData<SD> > vals;
        int ret;
        if ((ret = database.get(key, vals)) != 0)
        {
          std::cerr << __func__ << ": key not found\n" << k;
          exit(1);
        }

        eval.resize(vals.size());
        eval.resizeObs(vals[0].data().numElem());

        for(int i=0; i < vals.size(); ++i)
        {
          SD sval = vals[i].data();
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


  template<typename K, typename D>
    std::vector<D> RadmatDBInterface<K,D>::getEnsem(const std::vector<K> &kys) const
    {

      typedef typename ENSEM::EnsemScalar<D>::Type_t SD;

      POW2_ASSERT(m_haveparams);

      FILEDB::AllConfStoreDB<ADATIO::SerialDBKey<K> , ADATIO::SerialDBData<SD> > database;
      if (database.open(m_params.dbFile, O_RDONLY, 0400) != 0)
      {
        std::cerr << __func__ << ": error opening dbase= " << m_params.dbFile << std::endl;
        exit(1);
      }

      std::vector<D> data;

      try
      {
        typename std::vector<K>::const_iterator k;

        for(k = kys.begin(); k != kys.end(); ++k)
        {
          D eval;

          ADATIO::SerialDBKey<K> key;
          key.key() = *k;

          std::vector< ADATIO::SerialDBData<SD> > vals;
          int ret;
          if ((ret = database.get(key, vals)) != 0)
          {
            std::cerr << __func__ << ": key not found\n" << *k;
            exit(1);
          }

          eval.resize(vals.size());
          eval.resizeObs(vals[0].data().numElem());

          for(int i=0; i < vals.size(); ++i)
          {
            SD sval = vals[i].data();
            ENSEM::pokeEnsem(eval, sval, i);
          }

          data.push_back(eval);
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

      return data;

    }





} // namespace radmat











#endif /* RADMAT_DATABASE_INTERFACE_H */
