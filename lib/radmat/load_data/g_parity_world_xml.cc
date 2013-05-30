/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : g_parity_world_xml.cc

 * Purpose :

 * Creation Date : 24-04-2013

 * Last Modified : Wed May  8 09:50:51 2013

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "g_parity_world.h"
#include "formfac/formfac_qsq.h"
#include "hadron/irrep_util.h"
#include <sstream>


using namespace std;

namespace radmat
{

  namespace gParityWorld
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


      // stringy mom
      std::string stringmom(const ADATXML::Array<int> &inp)
      {
        ADATXML::Array<int> can = FF::canonicalOrder(inp);
        std::stringstream ss;
        ss << can[0] << can[1] << can[2];
        return ss.str(); 
      }


      // get the star
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


    } // anonomyous 


    // GParityContinuumStateXML
    //! write to a string for error
    std::string toString(const GParityContinuumStateXML &op)
    {
      std::stringstream ss;
      ss << "J = " << op.J << " H = ";
      for(int i =0; i < op.H.size(); ++i)
        ss << op.H[i] << " ";
      ss << " parity " << op.parity;
      ss << " mom = ";
      for(int i = 0; i < op.mom.size(); ++i)
        ss << "[" <<  stringmom(op.mom[i]) << "] ";
      ss << " twoI_z " << op.twoI_z << " op_stem " << op.op_stem << " create "; 
      ss << op.creation_op << " smear " << op.smearedP;
      return ss.str();
    }

    //! stream this thing
    std::ostream& operator<<(std::ostream &o , const GParityContinuumStateXML &e)
    {
      o << toString(e); 
      return o;
    }

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityContinuumStateXML &op)
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
      doXMLRead(ptop,"isProjected",op.isProjected,__PRETTY_FUNCTION__); 

      // update for empty helicities here.. if you dont want to consider one then explicitly provide a list
      if (op.H.size() == 0)
      {
        op.H.resize(2*op.J +1); 
        for(int h = -op.J; h < op.J +1; ++h)
          op.H[h + op.J] = h;
      }

      // fill the star if asked to
      if(op.fill_star)
        op.mom = star_p(op.mom);

    }

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityContinuumStateXML &op)
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
      write(xml,"isProjected",op.isProjected); 
      ADATXML::pop(xml);
    }



    // GParityListElement
    // stream
    std::string toString(const GParityInsertionXML::GParityListElement &e)
    {
      std::stringstream ss;
      ss << e.charge_coefficient << " X " << e.op_stem ;
      return ss.str(); 
    }

    std::ostream& operator<<(std::ostream &o, const  GParityInsertionXML::GParityListElement &e)
    {
      o << toString(e); 
      return o; 
    }

    //!  xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityInsertionXML::GParityListElement &e)
    {
      ADATXML::XMLReader ptop(xml,path); 
      doXMLRead(ptop,"op_stem",e.op_stem,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"charge_coefficient",e.charge_coefficient,__PRETTY_FUNCTION__);
    }

    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityInsertionXML::GParityListElement &e)
    { 
      ADATXML::push(xml,path);
      write(xml,"op_stem",e.op_stem);
      write(xml,"charge_coefficient",e.charge_coefficient);
      ADATXML::pop(xml);
    }

    // Insertion
    std::string toString(const GParityInsertionXML::Insertion &op)
    {
      std::stringstream ss;
      ss << "active = " << op.active << " J = " << op.J << " H = {";
      for(int i = 0; i < op.H.size(); ++i)
        ss << op.H[i] << " ";
      ss << "} parity = " << op.parity << " twoI_z = " << op.twoI_z
        << " create = " << op.creation_op
        << " smeared = " << op.smearedP << "\nphotons";
        for(int i = 0; i < op.photons.size(); ++i)
          ss << op.photons[i] << "\n";
      return ss.str();
    }

    std::ostream& operator<<(std::ostream &o, const  GParityInsertionXML::Insertion &e)
    {
      o << toString(e);
      return o; 
    }
    //!  xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityInsertionXML::Insertion &op)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"active",op.active,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"J",op.J,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"H",op.H,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"parity",op.parity,__PRETTY_FUNCTION__);      
      doXMLRead(ptop,"twoI_z",op.twoI_z,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"creation_op",op.creation_op,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"smearedP",op.smearedP,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"photons",op.photons,__PRETTY_FUNCTION__); 

      // update for empty helicities here.. if you dont want to consider one then explicitly provide a list
      if (op.H.size() == 0)
      {
        op.H.resize(2*op.J +1); 
        for(int h = -op.J; h < op.J +1; ++h)
          op.H[h + op.J] = h;
      }

    }

    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityInsertionXML::Insertion &op)
    {
      ADATXML::push(xml,path);
      write(xml,"active",op.active); 
      write(xml,"J",op.J);
      write(xml,"H",op.H);
      write(xml,"parity",op.parity);
      write(xml,"twoI_z",op.twoI_z);
      write(xml,"creation_op",op.creation_op);
      write(xml,"smearedP",op.smearedP);
      write(xml,"photons",op.photons);
      ADATXML::pop(xml);
    }

    // GParityInsertionXML
    std::string toString(const GParityInsertionXML &e)
    {
      std::stringstream ss; 
      ss << "time =" << e.time << "\nspace = " << e.space
        << " pmax = " << e.pmax << " t_slice = " << e.t_slice;
      return ss.str(); 
    }
    std::ostream& operator<<(std::ostream &o, const  GParityInsertionXML &e)
    {
      o << toString(e); 
      return o;
    }

    //!  xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityInsertionXML &e)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"pmax",e.pmax,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"pmin",e.pmin,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"t_slice",e.t_slice,__PRETTY_FUNCTION__); 
      doXMLRead(ptop,"time",e.time,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"space",e.space,__PRETTY_FUNCTION__); 
    }

    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityInsertionXML &e)
    {
      ADATXML::push(xml,path);
      write(xml,"pmax",e.pmax);
      write(xml,"pmin",e.pmin);
      write(xml,"t_slice",e.t_slice);
      write(xml,"time",e.time);
      write(xml,"space",e.space); 
      ADATXML::pop(xml);
    }

    //! write state to a string
    std::string toString(const GParityContinuumMatElemXML::State &e)
    {
      std::stringstream ss;
      ss << "state = " << e.state << " t_slice = " << e.t_slice;
      return ss.str();
    }

    //! stream a state
    std::ostream& operator<<(std::ostream &o , const GParityContinuumMatElemXML::State &e)
    {
      o << toString(e);
      return o; 
    }

    //! state xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityContinuumMatElemXML::State &e)
    {
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"state",e.state,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"t_slice",e.t_slice,__PRETTY_FUNCTION__); 
    }

    //! state xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityContinuumMatElemXML::State &e)
    {
      ADATXML::push(xml,path);
      write(xml,"state",e.state);
      write(xml,"t_slice",e.t_slice);
      ADATXML::pop(xml);
    }

    //! write it to a string
    std::string toString(const GParityContinuumMatElemXML &e)
    {
      std::stringstream ss;
      ss << "sink = " << e.sink;
      ss<< "\ninsertion = " << e.insertion;
      ss << "\nsource = " << e.source;
      return ss.str();
    }

    //! stream it
    std::ostream& operator<<(std::ostream &o, const GParityContinuumMatElemXML &e)
    {
      o << toString(e);
      return o;
    }

    //! xml reader
    void read(ADATXML::XMLReader &xml, const std::string &path, GParityContinuumMatElemXML &e)
    { 
      ADATXML::XMLReader ptop(xml,path);
      doXMLRead(ptop,"sink",e.sink,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"insertion",e.insertion,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"source",e.source,__PRETTY_FUNCTION__);
      doXMLRead(ptop,"ensemble",e.ensemble,__PRETTY_FUNCTION__); 
    }

    //! xml writer
    void write(ADATXML::XMLWriter &xml, const std::string &path, const GParityContinuumMatElemXML &e)
    { 
      ADATXML::push(xml,path);
      write(xml,"sink",e.sink);
      write(xml,"insertion",e.insertion);
      write(xml,"source",e.source);
      write(xml,"ensemble",e.ensemble);
      ADATXML::pop(xml);
    }

  } // gParityWorld


} // radmat
