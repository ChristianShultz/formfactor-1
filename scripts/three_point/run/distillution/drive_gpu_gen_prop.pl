#!/usr/bin/perl
#
# this file is suppose to be smart enough so I don't have to work hard
#  
#  at the time of writing it runs harom/chroma to get solution vectors/perambulators
#

use strict; 
use Params; 
use File::Basename;

my $hack = undef; 
# my $hack = "true"; 

die "usage: $0 <seqno> <delta_t> <output_path> " unless $#ARGV == 2; 

my $seqno = $ARGV[0];
my $outpath = $ARGV[2];
my $dt = $ARGV[1];

# an instance of the Param "class" to hold lattice specific stuff
my $param = Params->new(); 

# specify things that I want 
$param->seqno($seqno);
$param->num_vecs(64);       # max 
$param->num_space_dils(64); # effectively turn off distillution  
$param->num_time_dils(0);   # this variable is stupid
$param->nt_corr(40);        
$param->L_t(128);
$param->L_s(16);

# t_sources specifies sources
# delta_t specifies sinks 
#   this multiplies automatically and gets big quick 
$param->t_sources([0]);
$param->delta_t([$dt]);

# configuration stuff
$param->config_type("SCIDAC");
$param->stem("szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per");
$param->u_mass(-0.0743); 
$param->d_mass(-0.0743); 
$param->s_mass(-0.0743); 

# inversion stuff
$param->prop_dist_inverter("BICG_STAB");
$param->conn_quark_line_pool([1,2,3,4]);
$param->annih_quark_line_pool([1,2]);

# THIS IS WHAT SITS IN THE CFG.lime header and I assume it to be correct
$param->link_smear_type("STOUT_SMEAR");
$param->link_smear_fac(0.15625);
$param->link_smear_num(16);
$param->clover_c_s(1.58954576550983);
$param->clover_c_t(0.902866347885228);
$param->xi(3.5); # xi_mom also?
$param->xi_0(4.3);
$param->nu(1.265);
$param->rho(0.14);
$param->n_smear(2); 


# run everything from scratch
$param->scratchy("true");

# test lattice hackery 
if ( $hack ) 
{
  $param->skeleton("true");
  $param->num_vecs(4); 
  $param->num_space_dils(4);
  $param->nt_corr(16); 
  $param->L_t(16); 
  $param->L_s(4); 
  $param->delta_t([6]); 
  $param->config_type("WEAK_FIELD"); 
  $param->stem("fred"); 
  $param->prop_dist_inverter("CG_TEST"); 
  $param->link_smear_num(2); 
}


# be really noisy 
$param->verbose("true");

# gen a bunch of xml and get some file names
my $srcref =  $param->write_dist_src_ini_xml();
my $propref = $param->write_prop_dist_chroma_ini_xml(); 
my $peramref = $param->write_peram_dist_harom_ini_xml();

# deref some return values  , first value is the file name of the ini.xml , second is 
# a ref to an array of all the .mod files that we made assuming successfull completion 
# of the run code

my ($src,$srcs) = @{$srcref};
my ($prop,$props) = @{$propref};
my ($peram,$perams) = @{$peramref};

# try to run harom/chroma using the generated xml
#    first make the distillution source/noise V * noise
#    second get the solution vectors M^-1 * V * noise
#    third make some distillution perambulators noise' * V' * M^-1 * V * noise



&run_harom_ptx($src, $outpath."/".basename($src).".out" , " "); 
&run_chroma_ptx($prop, $outpath."/".basename($prop).".out" , " "); 
&run_harom_ptx($peram, $outpath."/".basename($peram).".out" , " "); 


# shove all the copy back junk into an array
# my @rcps = (); 
# push @rcps , @{$srcs} , @{$props} , @{$perams};
#
# copy lots of junk back to cache
# print "COPY BACK FILES\n";
# foreach my $ele (@rcps)
# {
#   print $ele . "\n";
# }
#
# $param->copy_back_scratch_rcp_arr(\@rcps);

if ( $param->scratchy() )
{
  my @sprops = (); 
  my @sperams = (); 

  push @sprops , @{$srcs} , @{$props};
  push @sperams , @{$perams};

  foreach my $el ( @sprops ) 
  {
    print "prop $el \n";
  }

  foreach my $el ( @sperams ) 
  {
    print "peram $el \n";
  }

  $param->scratchy_perams($param->scratchify_file_name_arr(\@sperams)); 
  $param->scratchy_props($param->scratchify_file_name_arr(\@sprops)); 
}

my $redstar = $param->write_redstar_xml();
my $harom_hadron_node = $param->write_harom_hadron_node_xml(); 


# Now try to run redstar_gen_graph
#
# Then run harom to generate the nodes
#
# Then run redstar_npt to multiply lots of stuff

&run_redstar_gen_graph($redstar,$outpath."/".basename($redstar).".gen");
&run_harom_ptx($harom_hadron_node,$outpath."/".basename($harom_hadron_node).".out","");
&run_redstar_npt($redstar,$outpath."/".basename($redstar).".npt");


# Now copy back npt_output.sdb, graph.sdb , and hadron_node.sdb

my $graphdb = $param->get_graph_db_string(); 
my $nptdb = $param->get_output_db_string(); 
my $hadnoddb = $param->get_hadron_node_output_db(); 

my $location = "/u/home/shultz/gen_props/szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per/zero_graph/${dt}/";

my @rcps = ();
my @perambulators = @{ $peram->scratchy_props() };
push @rcps , $graphdb , $nptdb , $hadnoddb, @perambulators;
foreach my $d (@rcps)
{
  my $base = basename($d); 
  my $db = $location . $base; 
  $param->copy_back_scratch_rcp($db);  
}


# clean up ?

# code -i inp -o out
#
#    delete inp and out 
#
unlink ( $outpath."/".$src.".out" ); 
unlink ( $outpath."/".$prop.".out" );
unlink ( $outpath."/".$peram.".out" );

unlink ( $src ); 
unlink ( $prop ); 
unlink ( $peram );




print "YOU WIN!\n";

exit ( 0 ) ; 




#
#
#   RUN CODE
#
#


sub mpi_info
{
  my $num_thread = `grep '^processor' /proc/cpuinfo | wc -l ` ;
  chomp $num_thread; 
  $num_thread /= 2; 

  my $nodefile = $ENV{'PBS_NODEFILE'};

  my $num_nodes = `cat $nodefile | wc -l` ;
  chomp $num_nodes; 

  my $base = dirname($0);
  my $numa = "${base}/numa_script.sh";

  my $mpi = $ENV{'MPIHOME'};

  my $run = "${mpi}/bin/mpirun_rsh -rsh -np $num_nodes -hostfile $nodefile MV2_ENABLE_AFFINITY=0 OMP_NUM_THREADS=2 $numa";


  print " MPI RUN INFO \n";
  print " num thread = $num_thread \n ";
  print " nodefile = $nodefile  \n" ;
  print " num_nodes = $num_nodes  \n" ;
  print " numa = $numa  \n" ;
  print " mpihome = $mpi \n";
  print " run = $run \n"; 


  die ( " cant find numa script : $numa " ) unless -e $numa;

  my @r = (); 
  push @r , $run , $num_nodes; 
  return \@r ; 
}


sub run_chroma_ptx
{
  my ($input_file , $output_file, $chroma_opts) = @_;

  my ($run,$num_nodes) = @{ &mpi_info() }; 

  $chroma_opts .= " -iogeom 1 1 1 $num_nodes  -geom 1 1 1 $num_nodes"; 

  my $chroma_exe = "/home/shultz/git-builds/chroma-jitptx/parscalar-Nd4/bin/chroma";
  my $chroma_command = "$run $chroma_exe $chroma_opts -i $input_file -o $output_file";
 
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($chroma_exe,"run_chroma_ptx");


  print " \n\n RUNNING CHROMA WITH \n\n $chroma_command \n\n";

  system ( " $chroma_command " ) == 0 || die (" some error for chroma, output : $output_file ");
}


sub run_harom_ptx
{
  my ($input_file , $output_file, $chroma_opts) = @_; 
  my $chroma_exe = "/home/shultz/git-builds/harom-jitptx/gpu-parscalar-Nd3/bin/harom";

  my ($run,$num_nodes) = @{ &mpi_info() }; 

  $chroma_opts .= " -iogeom 1 1 $num_nodes "; 

  my $chroma_command = "$run $chroma_exe $chroma_opts -i $input_file -o $output_file";
 
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($chroma_exe,"run_harom_ptx");

  print " \n\n RUNNING HAROM WITH \n\n $chroma_command \n\n";

  system ( " $chroma_command " ) == 0 || die (" some error for harom, output : $output_file \n $_ ");
}

sub run_redstar_npt
{
  my ($input_file , $output_file) = @_;
  my $exe = "/home/shultz/git-builds/redstar/bin/redstar_npt";
  my $cmd = "$exe $input_file $output_file ";
  
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($exe,"run_redstar_npt");

  system ( " $cmd " ) == 0 || die (" some error for redstar_npt, output : $output_file " ); 
}

sub run_redstar_gen_graph
{

  my ($input_file , $output_file, $options) = @_;
  my $exe = "/home/shultz/git-builds/redstar/bin/redstar_gen_graph";
  my $cmd = "$exe $input_file $output_file ";
  
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($exe,"run_redstar_gen_graph");

  system ( " $cmd " ) == 0 || die (" some error for redstar_npt, output : $output_file " ); 
}


sub check_exe
{
  my ($exe,$caller) = @_; 
  if ( ! -e $exe ) 
  {
    print " $exe does not exist! \n called by $caller \n";
    die("unholy death");
  }
}








