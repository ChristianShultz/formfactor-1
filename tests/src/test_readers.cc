/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 15-08-2012

 * Last Modified : Tue Oct 23 11:41:00 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include "radmat/load_data/simple_world.h"
#include "io/adat_xmlio.h"

namespace radmat
{

  //void make_test_file(const std::string &filename);

  tester test_readers(void)
  {
    tester m_test(__func__);

    simpleWorld::ContinuumStateXML csxml;

    csxml.J = 0;
    csxml.parity = true;
    csxml.mom.resize(1); 
    csxml.mom[0].resize(3);
    csxml.mom[0][0] = 1;
    csxml.mom[0][1] = 0;
    csxml.mom[0][2] = -1;
    csxml.twoI_z = 2;
    csxml.op_stem = "pion_pionxD0_J0__J0_";
    csxml.creation_op = true;
    csxml.smearedP = true;


    simpleWorld::ContinuumLorentzMatElem CLME, CLME_in;
    simpleWorld::ContinuumLorentzMatElem::State source,sink,ins;
    
    source.state = csxml;
    source.t_slice = 0;
    csxml.creation_op = false;
    sink.t_slice = 24;
    sink.state = csxml;
    csxml.op_stem = "b_b0xD0_J0__J0_";
    csxml.twoI_z = 0;
    ins.state = csxml;
    ins.t_slice = -3;

    CLME.source = source;
    CLME.sink = sink;
    CLME.lorentz.resize(4);
    CLME.lorentz[0] = ins;

    ins.state.op_stem = "rho_rhoxD0_J0__J1_";
    CLME.lorentz[1] = ins;
    CLME.lorentz[2] = ins;
    CLME.lorentz[3] = ins;


    ADATXML::XMLBufferWriter xml;
    simpleWorld::write(xml,"Key",CLME);
    
   // xml.print(std::cout);
    std::string xml_name("simpleWorldXMLout.xml");
    std::ofstream out(xml_name.c_str());
    xml.print(out);
    out.close();

    ADATXML::XMLReader xml_in(xml_name);
    ADATXML::XMLBufferWriter xml2;
    simpleWorld::read(xml_in,"/Key",CLME_in);
    simpleWorld::write(xml2,"Key",CLME_in);
    std::string xml_name2 = xml_name + std::string("_2");
    std::ofstream out2(xml_name2.c_str());
    xml2.print(out2);
    out2.close();

    // user should run diff on two output files if they care
    TESTER_TEST(m_test,true,"foobar");

    return m_test;
  }


} // namespace radmat


