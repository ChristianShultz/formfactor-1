#ifndef DISCONNECTED_GRAPH_NUKER_H
#define DISCONNECTED_GRAPH_NUKER_H


#include "radmat/utils/splash.h"
#include "radmat/driver/radmat_driver_props.h"

#include "io/adat_xmlio.h"

#include "adat/adat_stopwatch.h"

#include "hadron/ensem_filenames.h"
#include "hadron/hadron_npart_npt_conn_graph.h"

#include <vector>
#include <complex>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sstream>
#include <exception>



namespace radmat
{

  struct DisconnectedGraphNuker
  {

    void find_nukes(std::vector<Hadron::KeyHadronNPartNPtCorr_t> &keys,
        const std::string &hadron_node_xml );


    void dump_nukes(const std::string &filename);

    std::vector<Hadron::KeyHadronNPartNPtConnGraph_t> nuke_list;   

  };




} // namespace radmat





#endif /* DISCONNECTED_GRAPH_NUKER_H */
