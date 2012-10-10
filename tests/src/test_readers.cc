/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name :

 * Purpose :

 * Creation Date : 15-08-2012

 * Last Modified : Tue Sep 25 09:35:05 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/

#include "../headers/tester.h"
#include "../headers/test_common.h"
#include <stdio.h>
#include <fstream>
#include "radmat/load_data/read_redstar_key_list_xml.h"

namespace radmat
{

  //void make_test_file(const std::string &filename);

  tester test_readers(void)
  {
    tester m_test(__func__);

/*    std::string filename("junk.xml");

    make_test_file(filename);

    
    redstarKeyListXMLProp_t foobar = readKeys(filename);

    TESTER_TEST( m_test, remove(filename.c_str()) == 0, "something bad happend");
*/
    return m_test;
  }


/*
  void make_test_file(const std::string &filename)
  {

    std::ofstream out(filename.c_str());
  
//out << " <?xml version=\"1.0\"?> \n";
 out << " <Keys> <redstarKeyList> <elem> <NPoint> <elem> <t_slice>24</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> <elem> <t_slice>-3</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>rho_rhoxD0_J0__J1_T1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>3</row> <twoI_z>0</twoI_z> </Irrep> </elem> <elem> <t_slice>0</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>true</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> </NPoint><ensemble>fred</ensemble>";

 out << " </elem> <elem> <NPoint> <elem> <t_slice>24</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> <elem> <t_slice>-3</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>b_b0xD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>0</twoI_z> </Irrep> </elem> <elem> <t_slice>0</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>true</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> </NPoint> <ensemble>fred</ensemble> </elem>";

 out <<"<elem> <NPoint> <elem> <t_slice>24</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> <elem> <t_slice>-3</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>rho_rhoxD0_J0__J1_T1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>0</twoI_z> </Irrep> </elem> <elem> <t_slice>0</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>true</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> </NPoint> <ensemble>fred</ensemble> ";

out << " </elem> <elem> <NPoint> <elem> <t_slice>24</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> <elem> <t_slice>-3</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>rho_rhoxD0_J0__J1_T1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>false</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>2</row> <twoI_z>0</twoI_z> </Irrep> </elem> <elem> <t_slice>0</t_slice> <Irrep> <CGs> </CGs> <Operators> <elem> <name>pion_pionxD0_J0__J0_A1</name> <smear></smear> <mom_type>0 0 0</mom_type> </elem> </Operators> <creation_op>true</creation_op> <smearedP>true</smearedP> <mom>0 0 0</mom> <row>1</row> <twoI_z>2</twoI_z> </Irrep> </elem> </NPoint> <ensemble>fred</ensemble> </elem> </redstarKeyList> </Keys>";

    out.close();

  }

*/


} // namespace radmat


