#!/usr/bin/perl
#
#

use strict; 
use Params; 



my $param = Params->new(); 
$param->skeleton("true");
$param->seqno("10");
$param->npt_list_xml("npt.list.xml");
$param->stem("weak_field");
$param->num_vecs(10);
$param->nt_corr(16);
$param->t_sources([0]);
$param->L_t(16);
$param->L_s(4);
$param->num_space_dils(10); 
$param->num_time_dils(0);
$param->config_type("WEAK_FIELD");
$param->prop_dist_inverter("CG_TEST");
$param->conn_quark_line_pool([1,2,3,4]);
$param->annih_quark_line_pool([1,2]);


$param->write_dist_src_ini_xml();
$param->write_prop_dist_chroma_ini_xml(); 
$param->write_peram_dist_harom_ini_xml();
$param->write_redstar_xml();
$param->write_harom_hadron_node_xml(); 
