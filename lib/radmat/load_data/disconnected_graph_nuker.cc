/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : disconnected_graph_nuker.cc

 * Purpose :

 * Creation Date : 03-06-2013

 * Last Modified : Tue Jun  4 11:05:43 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/



#include "disconnected_graph_nuker.h"
#include "hadron/ensem_filenames.h"
#include "hadron/hadron_node.h"
#include "hadron/hadron_npart_npt_conn_graph.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "hadron/hadron_npart_irrep.h"
#include "hadron/hadron_timeslice.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"
#include <string>
#include <vector>
#include <map>
#include <algorithm>



/* so the npoint from the xml will have the insertion operator name
   the npoint from the grap db will also contain an operator name

   in the isoscalar case we get annihilation diagrams so we have to make sure
   the thing that we are zeroing actually came from the insertion.  This is 
   checked 2 ways.  

   1) insist that the orig_npt field in the graph key is the one we looked for
   2) insist that the name of the operator in the graph that we are killing 
   corresponds to one in the list of insertion names
 */



namespace radmat
{


  namespace 
  {


    // this is useful 
    template<typename K, typename V> 
      struct quick_map
      {

        void insert(const K &k , const V &v)
        {
          typename std::map<K,V>::const_iterator it; 
          it = data.find(k); 
          if( it == data.end() )
            data[k] = v; 
        }


        typename std::vector<V> vals(void) const
        {
          typename std::map<K,V>::const_iterator it; 
          typename std::vector<V> r;   
          for(it = data.begin(); it != data.end(); ++it)
            r.push_back(it->second); 
          return r; 
        }

        typename std::map<K,V> data;
      };


    // pick the irrep from an npoint 
    std::vector<Hadron::KeyHadronNPartIrrep_t> 
      unique_elem(const std::vector<Hadron::KeyHadronNPartNPtCorr_t> &keys,
          const int npt_elem)
      {
        quick_map<std::string , Hadron::KeyHadronNPartIrrep_t> mapp; 
        std::vector<Hadron::KeyHadronNPartNPtCorr_t>::const_iterator it; 
        for(it = keys.begin(); it != keys.end(); ++it)
          mapp.insert(Hadron::ensemFileName(it->npoint[npt_elem].irrep),
              it->npoint[npt_elem].irrep); 

        return mapp.vals(); 
      }


    // being lazy
    typedef Hadron::KeyHadronNPartNPtConnGraph_t GK; 
    typedef Hadron::ValTimeSlice_t GV;


    // open up the graph_db and get all of the keys
    std::vector<GK> graph_keys(const std::string &graph_db)
    {

      // Open DB
      FILEDB::ConfDataStoreDB< ADATIO::SerialDBKey<GK>, ADATIO::SerialDBData<GV> > database;
      if ( database.open(graph_db, O_RDONLY, 0400) != 0)
      {
        std::cerr << __PRETTY_FUNCTION__ << ": error opening " << graph_db << std::endl;
        exit(1); 
      }

      std::vector< ADATIO::SerialDBKey<GK> > keys;
      std::vector< ADATIO::SerialDBKey<GK> >::const_iterator it; 
      database.keys(keys);

      std::vector<GK> kys;

      for(it = keys.begin(); it != keys.end(); ++it)
        kys.push_back(it->key()); 

      return kys; 
    }


    // get the set of keys corresponing to an npoint function (n = 1 for me)
    std::vector<GK> find_npts(const int n, const std::string &graph_db)
    {
      std::vector<GK> all_keys = graph_keys(graph_db); 
      std::vector<GK> npt_keys;
      std::vector<GK>::const_iterator it; 


      for(it = all_keys.begin(); it != all_keys.end(); ++it)
        if(it->npoint.size() == n)
          npt_keys.push_back(*it); 

      return npt_keys; 
    }



    // unique names
    std::list<std::string> name_list(const std::vector<Hadron::KeyHadronNPartIrrep_t> &irreps)
    {
      std::list<std::string> names; 

      std::vector<Hadron::KeyHadronNPartIrrep_t>::const_iterator it; 
      for(it = irreps.begin(); it != irreps.end(); ++it)
        names.push_back(it->ops[1].name); // note stupid one particle op hardwire
      names.sort(); 
      names.unique(); 
      return names; 
    }


    // lazy
    std::list<std::string> name_list(const std::vector<Hadron::KeyHadronNPartNPtCorr_t> &v, 
        const int npt)
    {
      return name_list(unique_elem(v,npt)); 
    }


    // return true if we want to nuke this guy
    bool nuke_graph(const Hadron::KeyHadronNPartNPtConnGraph_t &g, 
        const int had_num,
        const std::list<std::string> &names)
    {
      if(g.npoint[1].orig_npt != had_num)
        return false; 

      std::list<std::string>::const_iterator it;
      it = std::find(names.begin(),names.end(),g.npoint[1].ops[1].name);

      return  !!!( it == names.end() );
    }


    // flip off smearing on the nuke list
    void turn_smearing_off(std::vector<Hadron::KeyHadronNPartNPtConnGraph_t> &v)
    {
      std::vector<Hadron::KeyHadronNPartNPtConnGraph_t>::iterator it; 
      for(it = v.begin(); it != v.end(); ++it)
        it->npoint[1].smearedP = false; 
    }


  } // namespace


  /*
     NOTE THE SERIES OF HARDWIRES FOR RADMAT HERE
   */




  void DisconnectedGraphNuker::find_nukes(std::vector<Hadron::KeyHadronNPartNPtCorr_t> &keys,
      const std::string &graph_db )
  {

    std::vector<Hadron::KeyHadronNPartNPtConnGraph_t> all_one_points = find_npts( 1 , graph_db ); 
    std::list<std::string> photon_names = name_list( keys, 2 );

    nuke_list.clear(); 

    std::vector<Hadron::KeyHadronNPartNPtConnGraph_t>::const_iterator it; 
    for(it = all_one_points.begin(); it != all_one_points.end(); ++it)
      if(nuke_graph(*it,2,photon_names))
        nuke_list.push_back(*it); 

    turn_smearing_off(nuke_list); 
  }


  void DisconnectedGraphNuker::dump_nukes(const std::string &fname)
  {
    ADATXML::XMLBufferWriter nuke_writer;
    ADATXML::Array<Hadron::KeyHadronNPartNPtConnGraph_t> nukes;
    nukes.resize(nuke_list.size());
    for(int i = 0; i < nuke_list.size(); ++i)
      nukes[i] = nuke_list[i];

    write(nuke_writer,"Keys",nukes); 
    std::ofstream out(fname.c_str());
    nuke_writer.print(out); 
    out.close(); 
  }


} // namespace radmat 
