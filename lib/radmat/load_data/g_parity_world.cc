/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : g_parity_world.cc

 * Purpose :

 * Creation Date : 25-04-2013

 * Last Modified : Wed 02 Oct 2013 01:20:23 PM EDT

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "g_parity_world.h"
#include "radmat/utils/pow2assert.h"
#include "hadron/irrep_util.h"
#include "formfac/formfac_qsq.h"

namespace radmat
{
  namespace gParityWorld
  {

    namespace
    {

      std::string string_mom(const ADATXML::Array<int> &b)
      {   
        std::stringstream ss; 
        for(int i = 0; i < b.size(); ++i)
          ss << b[i] << " ";
        return ss.str(); 
      }

      std::string string_mom_no_space(const ADATXML::Array<int> &b)
      {
        std::stringstream ss; 
        for(int i = 0; i < b.size(); ++i)
          ss << b[i];
        return ss.str(); 
      }


    }


    //! write it to a string for error
    std::string toString(const GParityState &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = " << op.H << " parity " << op.parity;
      ss << " mom = " << string_mom(op.mom);
      ss << " twoI_z " << op.twoI_z << " name " << op.name << " create "; 
      ss << op.creation_op << " smear " << op.smearedP;
      return ss.str();
    }

    //! stream a ContinuumExprPrim
    std::ostream& operator<<(std::ostream &o , const GParityState &e)
    {
      o << toString(e);
      return o;
    }

    std::string stateFileName(const GParityState &s)
    {
      std::stringstream ss; 
      ss << s.name << "_p" << string_mom_no_space(FF::canonicalOrder(s.mom))
        << ",J" << s.J << ",H" << s.H 
        << ",p" << string_mom_no_space(s.mom)
        << ",Iz" << s.twoI_z;
      return ss.str(); 
    }

    //! stream
    std::string toString(const GParityInsertion &e)
    {
      std::stringstream ss; 
      GParityInsertion::map_t::const_iterator it;
      for(it = e.begin(); it != e.end(); ++it)
        ss << it->first << " " << it->second << "\n";
      ss << "t_slice = " << e.t_slice << " create = " << e.create();
      return ss.str(); 
    }

    std::ostream& operator<<(std::ostream &o , const GParityInsertion &e)
    {
      o << toString(e);
      return o; 
    }


    //! stream 
    std::string toString(const GParityHelicityMatrixElement::State &e)
    {
      std::stringstream ss; 
      ss << "state = " << e.state << " t_slice = " << e.t_slice; 
      return ss.str(); 
    }

    std::ostream& operator<<(std::ostream &o, const GParityHelicityMatrixElement::State &e)
    {
      o << toString(e);
      return o; 
    }

    std::string stateFileName(const GParityHelicityMatrixElement::State &e)
    {
      std::stringstream ss; 
      ss << stateFileName(e.state) << ",t" << e.t_slice;
      return ss.str(); 
    }

    // GParityHelicityMatrixElement
    std::string toString(const GParityHelicityMatrixElement &e)
    {
      std::stringstream ss;
      ss << "sink = " << e.sink; 
      ss << "\ninsertion = " << e.insertion;
      ss << "\nsource = " << e.source;
      ss << "\nensemble = " << e.ensemble;
      return ss.str(); 
    }
    std::ostream& operator<<(std::ostream &o , const GParityHelicityMatrixElement &e)
    {
      o << toString(e);
      return o;
    }


    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////


    namespace
    {

      // make one of the internal states on the fly
      GParityState makeGParityState(const int J, const int H, const bool parity, const ADATXML::Array<int> &mom,
          const int twoI_z, const std::string name, const bool creation_op, const bool smearedP, const bool isProjected)
      {
        POW2_ASSERT(abs(H) <= J); // sanity
        GParityState ret; 
        ret.J = J;
        ret.H = H;
        ret.parity = parity;
        ret.mom = mom;
        ret.twoI_z = twoI_z;
        ret.name = name;
        ret.creation_op = creation_op;
        ret.smearedP = smearedP; 
        ret.isProjected = isProjected; 
        return ret; 
      }

      // loop over whatever the xml gave us to make a list of possible states
      std::vector<GParityState> getGParityStateFromXML(const GParityContinuumStateXML &xml)
      {
        std::vector<GParityState> ret; 
        for(int h = 0; h < xml.H.size(); ++h)
          for(int p = 0; p < xml.mom.size(); ++p)
            ret.push_back(makeGParityState(xml.J,xml.H[h],xml.parity,xml.mom[p],
                  xml.twoI_z,xml.op_stem,xml.creation_op,xml.smearedP,xml.isProjected));
        return ret; 
      }

      // photon is a list of coeff times a GParityState -- generate one
      GParityInsertion::photon getPhoton(const GParityInsertionXML::Insertion &ins, 
          const int H, const ADATXML::Array<int> &mom) 
      {
        GParityInsertion::photon ret; 

        for(int i = 0; i < ins.photons.size(); ++i)
          ret = ret + GParityInsertion::photon_frag(ins.photons[i].charge_coefficient,
              makeGParityState(ins.J,H,ins.parity,mom,ins.photons[i].twoI_z,ins.photons[i].op_stem,
                ins.creation_op,ins.smearedP,false) // insertions are not projected operators 
              );

        return ret; 
      }

      // sanity and get the answer
      bool getCreationOp(const GParityInsertionXML &e)
      {
        bool create; 

        // figure out if its a creation op
        if(e.time.active && !!!e.space.active)
          create = e.time.creation_op;
        else if (!!!e.time.active && e.space.active)
          create = e.space.creation_op;
        else if(e.time.creation_op == e.space.creation_op)
          create = e.time.creation_op;
        else
        {
          std::cerr << __func__ << "error selecting the creation op, exiting";
          exit(1);
        }

        return create; 
      }


      // determine the momentum transfer
      ADATXML::Array<int> momentumTransfer(const GParityState &sink, 
          const GParityState & source, 
          const bool create)
      { 
        ADATXML::Array<int> ret;
        ret.resize(3);

        ret[0] = source.mom[0] - sink.mom[0];
        ret[1] = source.mom[1] - sink.mom[1];
        ret[2] = source.mom[2] - sink.mom[2];

        if (create) 
          return -ret;

        return ret; 
      }


      // do the cut
      bool cutMomentum(const ADATXML::Array<int> &mom, const int minmom , const int maxmom)
      {
        int sq(0);
        for(int i = 0; i < 3; ++i)
          sq += mom[i]*mom[i];
        return   !!! ( (sq >= minmom) && (sq <= maxmom) ) ;  // false if its too big or too small 
      }


      bool acceptable_mom(const ADATXML::Array<int> &mom)
      {
        std::string lg = Hadron::generateLittleGroup(mom); 
        if(lg == "Oh")
          return true;
        if(lg == "D4")
          return true; 
        if(lg == "D2")
          return true; 
        if(lg == "D3")
          return true; 

        std::cout << __func__ << ": not supporting LG " << lg << std::endl;
        return false;
      }


      // generate the internal representation of the insertion from xml 
      std::pair<bool,GParityInsertion> makeInsertion(const GParityState &sink,
          const GParityInsertionXML &xml,
          const GParityState &source)
      {         
        GParityInsertion ret; 
        bool creation_op = getCreationOp(xml); 
        ADATXML::Array<int> mom = momentumTransfer(sink,source,creation_op);


        ret.set_create(creation_op); 
        ret.t_slice = xml.t_slice; 
        ret.mom = mom;

        // break early if we can't make it
        if(cutMomentum(mom,xml.pmin,xml.pmax))
          return std::pair<bool,GParityInsertion>(false,ret);

        if(!!!acceptable_mom(mom)) 
          return std::pair<bool,GParityInsertion>(false,ret);

        // do the temporal insertion
        if(xml.time.active) 
          ret.insert("t",getPhoton(xml.time,0,mom)); 

        // do the spatial insertion
        if(xml.space.active)
        {
          std::map<int,std::string> foo;
          foo[1] = "p";
          foo[0] = "0";
          foo[-1] = "m";

          for(int h = 0; h < xml.space.H.size(); ++h)
          {
            if(foo.find(xml.space.H[h]) != foo.end())
            {
              std::string lab = foo[xml.space.H[h]];
              ret.insert(lab,getPhoton(xml.space,xml.space.H[h],mom)); 
            }
            else
            {
              std::cerr << __func__ << ": error unknown helicity " << xml.space.H[h] << std::endl;
              exit(1); 
            }
          }
        }

        return std::pair<bool,GParityInsertion>(true,ret); 
      }


      std::vector<GParityHelicityMatrixElement::State> 
        getStates(const GParityContinuumMatElemXML::State &s)
        {
          std::vector<GParityState> e = getGParityStateFromXML(s.state);
          std::vector<GParityHelicityMatrixElement::State> ret; 
          std::vector<GParityState>::const_iterator it; 
          for(it = e.begin(); it != e.end(); ++it)
            ret.push_back(GParityHelicityMatrixElement::State(*it,s.t_slice)); 
          return ret; 
        }

      GParityHelicityMatrixElement makeMatrixElement(
          const GParityHelicityMatrixElement::State &sink, 
          const GParityInsertion &ins,
          const GParityHelicityMatrixElement::State &source,
          const std::string &ensemble)
      {
        GParityHelicityMatrixElement ret;
        ret.sink = sink;
        ret.insertion = ins;
        ret.source = source;
        ret.ensemble = ensemble; 
        return ret; 
      }

    } // anonymous 




    std::vector<GParityHelicityMatrixElement> 
      getGParityHelicityMatrixElementFromXML(const GParityContinuumMatElemXML &xml)
      {
        std::vector<GParityHelicityMatrixElement> ret; 
        std::vector<GParityHelicityMatrixElement::State> sinks,sources;  
        std::vector<GParityHelicityMatrixElement::State>::const_iterator sink, source; 
        sources = getStates(xml.source);
        sinks = getStates(xml.sink); 

        for(sink = sinks.begin(); sink != sinks.end(); ++sink)
          for(source = sources.begin(); source != sources.end(); ++source)
          {
            std::pair<bool,GParityInsertion> ins;
            ins = makeInsertion(sink->state, xml.insertion, source->state);
            if(ins.first)
              ret.push_back(makeMatrixElement(*sink,ins.second,*source,xml.ensemble)); 
          }

        if(ret.empty())
        {
          std::cerr << __func__ << ": error, returning an empty vector, aborting" << std::endl;
          exit(1);
        }

        return ret; 
      }



  } // gParityWorld

} // radmat
