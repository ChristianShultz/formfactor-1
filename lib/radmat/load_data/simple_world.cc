/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : simple_world.cc

 * Purpose :

 * Creation Date : 19-10-2012

 * Last Modified : Wed Jan 23 14:48:16 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "hadron/irrep_util.h"
#include "radmat/utils/pow2assert.h"
#include "formfac/formfac_qsq.h"
#include "simple_world.h"
#include <sstream>
#include <omp.h>


using namespace ENSEM;
using namespace ADATXML;
using namespace Hadron;
using namespace radmat::simpleWorld;

namespace radmat
{

  namespace simpleWorld
  {


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


        std::string stringmom(const ADATXML::Array<int> &inp)
        {
          ADATXML::Array<int> can = FF::canonicalOrder(inp);
          std::stringstream ss;
          ss << can[0] << can[1] << can[2];
          return ss.str(); 
        }


        ADATXML::Array<ADATXML::Array<int> > star_p(const ADATXML::Array<ADATXML::Array<int> > &pin)
        {

          std::list<ADATXML::Array<int> > work;
          std::map<std::string,std::list<ADATXML::Array<int> > > seen;
          std::map<std::string,std::list<ADATXML::Array<int> > >::const_iterator it;
          std::list<ADATXML::Array<int> >::const_iterator listit; 

          int nele(0); 

          const unsigned int sz = pin.size(); 

          for(unsigned int elem = 0; elem < sz; ++elem)
          {
            std::string p = stringmom(pin[elem]);
            it = seen.find(p);
            if (it == seen.end())
            {
              work = Hadron::generateLittleGroupMom(Hadron::generateLittleGroup(pin[elem]),
                                                    FF::canonicalOrder(pin[elem]));
              nele += work.size();
              seen.insert(std::pair<std::string,std::list<ADATXML::Array<int> > >(p,work)); 
            }        
          }


          ADATXML::Array<ADATXML::Array<int> > ret(nele);
          nele = 0;

          for(it = seen.begin(); it != seen.end(); ++it)
            for(listit = it->second.begin(); listit != it->second.end(); ++listit)
            {
              ret[nele] = *listit; 
              ++nele; 
            }

          return ret; 
        }


    } // namespace anonomyous 


    // ContinuumStatePrimitive
    //------------------------------------------------------------------------------------------

    //! write it to a string for error
    std::string toString(const ContinuumStatePrimitive &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = " << op.H << " parity " << op.parity;
      ss << " mom = " << op.mom[0] << " " << op.mom[1] << " " << op.mom[2];
      ss  << " twoI_z " << op.twoI_z << " name " << op.name << " create "; 
      ss << op.creation_op << " smear " << op.smearedP;
      return ss.str();
    }

    //! stream a ContinuumExprPrim
    std::ostream& operator<<(std::ostream& o , const ContinuumStatePrimitive &op)
    {
      o << toString(op);
      return o;
    }

  
    std::string stateFileName(const ContinuumStatePrimitive &s)
    {
      std::stringstream ss; 
      ss << s.name << ",J" << s.J << ",H" << s.H << ",p" << s.mom[0] << s.mom[1] << s.mom[2] << ",Iz" << s.twoI_z;
      return ss.str(); 
    }

    // ContinuumStateXML
    //------------------------------------------------------------------------------------------


    //! write to a string for error
    std::string toString(const ContinuumStateXML &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = {";
      for(int i = 0; i < op.H.size(); ++i)
        ss << op.H[i] << " ";
      ss << "} parity = " << op.parity << " mom = { ";
      ss << " fill_star " << op.fill_star;  

      for(int i = 0; i < op.mom.size(); ++i)
        ss <<  "(" <<op.mom[i][0] << " " << op.mom[i][1] << " " << op.mom[i][2] << ") ";
      ss << "} twoI_z = " << op.twoI_z << " op_stem " << op.op_stem << " create = "
        << op.creation_op << " smeared = " << op.smearedP;
      return ss.str();
    }

    //! stream this thing
    std::ostream& operator<<(std::ostream& o , const ContinuumStateXML &op)
    {
      o << toString(op);
      return o;
    }

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumStateXML &op)
    {   
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"J",op.J,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"H",op.H,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"parity",op.parity,__PRETTY_FUNCTION__);  
      doXMLRead(ptop,"fill_star",op.fill_star,__PRETTY_FUNCTION__);     
      doXMLRead(ptop,"mom",op.mom,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"twoI_z",op.twoI_z,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"op_stem",op.op_stem,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"creation_op",op.creation_op,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"smearedP",op.smearedP,__PRETTY_FUNCTION__);


      // update for empty helicities here.. if you dont want to consider one then explicitly provide a list
      if (op.H.size() == 0)
      {
        op.H.resize(2*op.J +1); 
        for(int h = -op.J; h < op.J +1; ++h)
          op.H[h + op.J] = h;
      }

      if(op.fill_star)
        op.mom = star_p(op.mom);
    } 

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumStateXML &op)
    {
      ADATXML::push(xml,path);
      write(xml,"J",op.J);
      write(xml,"H",op.H);
      write(xml,"parity",op.parity);
      write(xml,"mom",op.mom);
      write(xml,"twoI_z",op.twoI_z);
      write(xml,"op_stem",op.op_stem);
      write(xml,"creation_op",op.creation_op);
      write(xml,"smearedP",op.smearedP);
      ADATXML::pop(xml);
    }


    // ContinuumInsertion
    //------------------------------------------------------------------------------------------

    std::string toString(const ContinuumInsertion & ins) 
    {
      std::stringstream ss;
      std::map<std::string,ContinuumInsertion::op_insertion>::const_iterator it;
      ss << "t_slice = " << ins.t_slice << "\n";
      for(it = ins.insertion_map.begin(); it != ins.insertion_map.end(); ++it)
        ss << it->first << " = " <<  it->second << "\n";
      return ss.str(); 
    }

    std::ostream& operator<<(std::ostream & o , const ContinuumInsertion &ins)
    {
      o << toString(ins);
      return o; 
    }



    // ContinuumInsertionXML
    //------------------------------------------------------------------------------------------

    //! write to a string for error
    std::string toString(const ContinuumInsertionXML::Insertion &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = {";
      for(int i = 0; i < op.H.size(); ++i)
        ss << op.H[i] << " ";
      ss << "} parity = " << op.parity << " twoI_z = " << op.twoI_z
        << " op_stem " << op.op_stem << " create = " << op.creation_op
        << " smeared = " << op.smearedP;
      return ss.str();
    }

    //! stream this thing
    std::ostream& operator<<(std::ostream& o , const ContinuumInsertionXML::Insertion &op)
    {
      o << toString(op);
      return o;
    }

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumInsertionXML::Insertion &op)
    {   
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"J",op.J,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"H",op.H,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"parity",op.parity,__PRETTY_FUNCTION__);      
      doXMLRead(ptop,"twoI_z",op.twoI_z,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"op_stem",op.op_stem,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"creation_op",op.creation_op,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"smearedP",op.smearedP,__PRETTY_FUNCTION__);

      // update for empty helicities here.. if you dont want to consider one then explicitly provide a list
      if (op.H.size() == 0)
      {
        op.H.resize(2*op.J +1); 
        for(int h = -op.J; h < op.J +1; ++h)
          op.H[h + op.J] = h;
      }
    } 

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumInsertionXML::Insertion &op)
    {
      ADATXML::push(xml,path);
      write(xml,"J",op.J);
      write(xml,"H",op.H);
      write(xml,"parity",op.parity);
      write(xml,"twoI_z",op.twoI_z);
      write(xml,"op_stem",op.op_stem);
      write(xml,"creation_op",op.creation_op);
      write(xml,"smearedP",op.smearedP);
      ADATXML::pop(xml);
    }

    //! write to a string for error
    std::string toString(const ContinuumInsertionXML &op)
    {
      std::stringstream ss;
      ss << "t_slice = " << op.t_slice << "\ntime = " << op.time << "\nspace = " << op.space; 
      return ss.str();
    }

    //! stream this thing
    std::ostream& operator<<(std::ostream& o , const ContinuumInsertionXML &op)
    {
      o << toString(op); 
      return o;
    }

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumInsertionXML &op)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"pmax",op.pmax,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"t_slice",op.t_slice,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"time",op.time,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"space",op.space,__PRETTY_FUNCTION__);
    }

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumInsertionXML &op)
    {
      ADATXML::push(xml,path);
      write(xml,"t_slice",op.t_slice);
      write(xml,"time",op.time);
      write(xml,"space",op.space);
      ADATXML::pop(xml);
    }


    // ContinuumMatElem
    //------------------------------------------------------------------------------------------


    //! write state to a string
    std::string toString(const ContinuumMatElem::State &s)
    {
      std::stringstream ss;
      ss << "state = " << s.state << " t_slice = " << s.t_slice;
      return ss.str();
    }

    std::ostream& operator<<(std::ostream& o, const ContinuumMatElem::State &s)
    {
      o << toString(s);
      return o;
    }



    //! write it to a string for error
    std::string toString(const ContinuumMatElem & s)
    {
      std::stringstream ss;
      ss << " source = " << s.source << "\n";
      ss << "insertion = " << s.insertion << "\n";
      ss << "sink = " << s.sink << "\n";
      ss << "ensemble = " << s.ensemble; 
      return ss.str();
    }

    //! stream a ContinuumMatElem
    std::ostream& operator<<(std::ostream& o, const ContinuumMatElem & s)
    {
      o << toString(s);
      return o;
    }

    std::string stateFileName(const ContinuumMatElem::State &s)
    {
      std::stringstream ss;
      ss << stateFileName(s.state) << ",t" << s.t_slice;
      return ss.str(); 
    }

    // ContinuumMatElemXML
    //------------------------------------------------------------------------------------------



    //! write it to a string
    std::string toString(const ContinuumMatElemXML &s)
    {
      std::stringstream ss;
      ss << "source = " << s.source << "\n";
      ss << "insertion " << s.insertion << "\n";
      ss << "sink = " << s.sink << "\n";
      ss << "ensemble = " << s.ensemble; 
      return ss.str();
    }

    //! stream it
    std::ostream& operator<<(std::ostream& o, const ContinuumMatElemXML &s)
    {
      o << toString(s);
      return o;
    }


    //! write state to a string
    std::string toString(const ContinuumMatElemXML::State &s)
    {
      std::stringstream ss;
      ss << "state = " << s.state << " t_slice = " << s.t_slice;
      return ss.str(); 
    }

    std::ostream& operator<<(std::ostream& o , const ContinuumMatElemXML::State &s)
    {
      o << toString(s);
      return o;
    }


    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumMatElemXML::State &op)
    {   
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"t_slice",op.t_slice,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"state",op.state,__PRETTY_FUNCTION__);

    }

    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumMatElemXML::State &op)
    { 
      ADATXML::push(xml,path);
      write(xml,"state",op.state);
      write(xml,"t_slice",op.t_slice);
      ADATXML::pop(xml);
    }



    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, ContinuumMatElemXML &op)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"source",op.source,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"insertion",op.insertion,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"sink",op.sink,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"ensemble",op.ensemble,__PRETTY_FUNCTION__); 
    }

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const ContinuumMatElemXML &op)
    {
      ADATXML::push(xml,path);
      write(xml,"source",op.source);
      write(xml,"insertion",op.insertion);
      write(xml,"sink",op.sink);
      write(xml,"ensemble",op.ensemble); 
      ADATXML::pop(xml);
    }



    namespace
    {
      std::vector<ContinuumStatePrimitive> getContinuumStatePrimitiveFromXML(const ContinuumStateXML &xmlin)
      {
        ContinuumStatePrimitive base; 
        base.J = xmlin.J;
        base.parity = xmlin.parity;
        base.twoI_z = xmlin.twoI_z;
        base.name = xmlin.op_stem;
        base.creation_op = xmlin.creation_op;
        base.smearedP = xmlin.smearedP;
        base.mom.resize(3);


        std::vector<ContinuumStatePrimitive> ret; 


        for(int p = 0; p < xmlin.mom.size(); ++p)
        {
          POW2_ASSERT(xmlin.mom[p].size() == 3); 

          base.mom[0] = xmlin.mom[p][0];
          base.mom[1] = xmlin.mom[p][1];
          base.mom[2] = xmlin.mom[p][2];

          // allow for leaving the helicity bit blank so we dont always have to type them all in
          if(xmlin.H.size() == 0)
          {
            for(int h = -base.J; h < base.J + 1; ++h)
            {
              base.H = h;
              ret.push_back(base); 
            }
          }
          else
          {
            for(int h = 0; h < xmlin.H.size(); ++h)
            {
              base.H = xmlin.H[h];
              ret.push_back(base);
            }
          }

        } // loop momenta

#if 0
        std::vector<ContinuumStatePrimitive>::const_iterator it;
        std::cout << __func__ << std::endl;
        for(it = ret.begin(); it != ret.end(); ++it)
          std::cout << *it << std::endl;
#endif

        return ret; 
      }


      ContinuumStatePrimitive makeTemporalInsertion(const ContinuumInsertionXML::Insertion &time)
      { 
        POW2_ASSERT(time.J == 0);
        ContinuumStatePrimitive ret; 
        ret.J = 0;
        ret.H = 0;
        ret.parity = time.parity;
        ret.twoI_z = time.twoI_z;
        ret.name = time.op_stem;
        ret.creation_op = time.creation_op;
        ret.smearedP = time.smearedP;
        return ret; 
      } 

      namespace
      {
        struct triplet
        {
          triplet(const ContinuumInsertionXML::Insertion &space)
          {
            POW2_ASSERT(space.J == 1);
            plus.J = 1;
            plus.parity = space.parity;
            plus.twoI_z = space.twoI_z;
            plus.name = space.op_stem;
            plus.creation_op = space.creation_op;
            plus.smearedP = space.smearedP; 
            minus = plus;
            zero = plus; 
            plus.H = 1;
            minus.H = -1;
            zero.H = 0; 
          }
          ContinuumStatePrimitive plus,minus,zero;
        };
      } // namespace anonomyous 


      triplet makeCircularSpatialInsertion(const ContinuumInsertionXML::Insertion &space)
      {
        triplet circ(space); 
        return circ; 
      } 


      ContinuumInsertion getContinuumInsertionFromXML(const ContinuumInsertionXML &xmlin)
      {
        ContinuumInsertion ret;
        ret.t_slice = xmlin.t_slice; 
        ContinuumStatePrimitive t = makeTemporalInsertion(xmlin.time);
        triplet circ = makeCircularSpatialInsertion(xmlin.space); 

        ret.insert("t",t);


        // allow for leaving the helicity bit blank so we dont always have to type them all in
        ContinuumInsertionXML my_insertion(xmlin);

        if(my_insertion.space.H.size() == 0)
        {
          my_insertion.space.H.resize(2*my_insertion.space.J + 1);

          for(int h = -my_insertion.space.J ; h < my_insertion.space.J +1; ++h)
            my_insertion.space.H[h + my_insertion.space.J] = h;
        }

        // hardwire for J = 1 xml input type
        // allow for not considering certain bits of the insertion 
        for(int h = 0; h < xmlin.space.H.size(); ++h)
        {
          if(my_insertion.space.H[h] == -1)
          {
            ret.insert("m",circ.minus);
            continue;
          }
          else if (my_insertion.space.H[h] == 0)
          {
            ret.insert("0",circ.zero);
            continue;
          }
          else if (my_insertion.space.H[h] == 1)
          {
            ret.insert("p",circ.plus);
          }
        }

        return ret; 
      }


      std::vector<ContinuumMatElem::State> getMatElemState(const ContinuumMatElemXML::State &xmlin)
      {
        std::vector<ContinuumMatElem::State> ret; 
        std::vector<ContinuumStatePrimitive> states =  getContinuumStatePrimitiveFromXML(xmlin.state); 
        std::vector<ContinuumStatePrimitive>::const_iterator it; 

        for(it = states.begin(); it != states.end(); ++it)
          ret.push_back(ContinuumMatElem::State(*it,xmlin.t_slice)); 

        return ret;
      }


      ADATXML::Array<int> getMomentumTransfer(const ContinuumMatElem::State &source, const ContinuumMatElem::State &sink)
      {
        ADATXML::Array<int> ret;
        ret.resize(3);
        // NB: need that extra minus sign or else all of the correlators are zero b/c of phase convention chosen for the insertion
        // to be daggered (un-daggered.. can't remember which but this is the way it needs to be)
        ret[0] = -(source.state.mom[0] - sink.state.mom[0]);
        ret[1] = -(source.state.mom[1] - sink.state.mom[1]);
        ret[2] = -(source.state.mom[2] - sink.state.mom[2]);
        return ret; 
      } 


    } // namespace anonomyous 



    std::vector<ContinuumMatElem> getContinuumMatElemFromXML(const ContinuumMatElemXML &xmlin)
    {
      std::vector<ContinuumMatElem> ret; 
      std::vector<ContinuumMatElem::State> source,sink;
      std::vector<ContinuumMatElem::State>::const_iterator it_source, it_sink;
      ContinuumInsertion ins; 

      source = getMatElemState(xmlin.source);
      sink = getMatElemState(xmlin.sink);
      ins = getContinuumInsertionFromXML(xmlin.insertion);

      for(it_source = source.begin(); it_source != source.end(); ++it_source)
        for(it_sink = sink.begin(); it_sink != sink.end(); ++it_sink)
        {
          ContinuumMatElem dum;
          dum.source = *it_source;
          dum.sink = *it_sink;
          dum.insertion = ins;
          dum.ensemble = xmlin.ensemble;

          ADATXML::Array<int> q = getMomentumTransfer(*it_source,*it_sink);
  
          int qq = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];
  
          // skip the guys for qq_space > pmax 
          if(qq > xmlin.insertion.pmax)
            continue; 

          std::map<std::string,ContinuumInsertion::op_insertion>::iterator it;

          for(it = dum.insertion.insertion_map.begin(); it != dum.insertion.insertion_map.end(); ++it)
            it->second.mom = q;

          ret.push_back(dum);
        }
      return ret; 
    }


  } // namespace simpleWorld

} // namespace radmat


