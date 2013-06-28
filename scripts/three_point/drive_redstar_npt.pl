#!/usr/bin/perl
#
#run run_redstar.csh


use strict; 
use Params; 
use File::Copy; 


die "usage: $0 <seqno>\n" unless @ARGV == 1; 

my $seqno = shift(@ARGV); 

my $npt_list = "npt.list.xml";
my $shell_script = "run_redstar.csh";


my $param = Params->new();
$param->seqno($seqno);
$param->npt_list($npt_list);

mkdir $seqno || die("Couldn't mkdir $seqno \nError: $_");

my $meson_hadron_node_file = $param->write_meson_hadron_node_xml();
my $redstar_xml_file = $param->write_redstar_xml(); 

move("$meson_hadron_node_file","${seqno}/${meson_hadron_node_file}") || die ("Couldn't mv $meson_hadron_node_file \nError: $_"); 
move("$redstar_xml_file","${seqno}/${redstar_xml_file}") || die("Couldn't mv $redstar_xml_file \n Error:$_"); 

chdir $seqno; 

# clean up anything old
my $had_graph = $param->hadron_node_graph();
my $had_node_out = $param->hadron_node_output_db();
my $out_db = $param->output_db(); 

if(-e $had_graph)
{
  unlink($had_graph);
}

if(-e $had_node_out)
{
  unlink($had_node_out);
}

if(-e $out_db)
{
  unlink($out_db); 
}

my $command = "$shell_script $redstar_xml_file"; 
print $command . "\n";

system ($command) == 0 || die($_); 

chdir "..";

return 0; 
