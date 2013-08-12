#!/usr/bin/perl
#
# this file is suppose to be smart enough so I don't have to work hard
#  
#  at the time of writing it runs harom/chroma to get solution vectors/perambulators
#

use strict; 
use ParamsDistillation; 
use File::Basename;

my $hack = undef; 
# my $hack = "true"; 

die "usage: $0 <seqno> <delta_t> <output_path> " unless $#ARGV == 2; 

my $seqno = $ARGV[0];
my $outpath = $ARGV[2];
my $dt = $ARGV[1];

# an instance of the Param "class" to hold lattice specific stuff
my $param = ParamsDistillation->new(); 

# REDSTAR AT VERSION 5 -- allow convertUDtoS
# VERSION 1 does not exist must be [2,5]
# $param->version(5); 
$param->redstar_version(6); 


# Distillation should be at version 1
$param->harom_version(1); 

# specify things that I want 
$param->seqno($seqno);
$param->num_vecs(64);       # max 
$param->nt_corr($dt + 1);  # this should avoid yet another way of mucking it up   
$param->L_t(128);         
$param->L_s(16);

# t_sources specifies sources
# delta_t specifies sinks 
#   this multiplies automatically and gets big quick 
$param->t_sources([0]);
$param->delta_t([$dt]);

# configuration stuff
$param->config_type("SZINQIO");
my $lattice = "szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per";
$param->stem($lattice);
my @proj_op = ("/u/home/shultz/optimized_operators/${lattice}/weights.${lattice}.xml");
$param->proj_op_xml(\@proj_op);
$param->convertUDtoL("true");   # this should be false if UD->S is true
$param->convertUDtoS("false");    # requires redstar version 5 or better 
$param->u_mass(-0.0743); 
$param->d_mass(-0.0743); 
$param->s_mass(-0.0743); 

# inversion stuff
$param->inverter("BICGSTAB");  # QUDA isn't working

# THIS IS WHAT SITS IN THE CFG.lime header and I assume it to be correct
$param->link_smear_type("STOUT_SMEAR");
$param->link_smear_fac(0.15625);
$param->link_smear_num(16);               # this is insane.. ask RGE
$param->clover_c_s(1.58954576550983);
$param->clover_c_t(0.902866347885228);
$param->xi(3.5); # xi_mom also?
$param->xi_0(4.3);
$param->nu(1.265);
$param->rho(0.14);
$param->n_smear(2); 


# run everything from scratch -- avoid any copying to  
#     lustre until we have a gen prop / perambulators
$param->scratchy(undef);

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
$param->diagnostic_level(1); 

#    
#   THIS IS THE OLD GEN PROP JUNK THAT WORKS WITH DISTILLATION
#
#   #    :M,N s/^/#  /  -- or switch to uncomment 
#   # gen a bunch of xml and get some file names
#   my $srcref =  $param->write_harom_distillation_source_ini_xml();
#   my $propref = $param->write_chroma_distillation_prop_ini_xml(); 
#   
#   
#   # deref some return values  , first value is the file name of the ini.xml , second is 
#   # a ref to an array of all the .mod files that we made assuming successfull completion 
#   # of the run code
#   
#   my ($src,$srcs) = @{$srcref};
#   my ($prop,$props) = @{$propref};
#   
#   # get some solution vectors up in /scratch
#   &run_harom_ptx($src, $outpath."/".basename($src).".out" , " "); 
#   &run_chroma_ptx($prop, $outpath."/".basename($prop).".out" , " "); 
#   
#   # need to do some redstar_gen_graph here to get the unsmeared list
#   my $redxml = $param->write_redstar_xml(); 
#   &run_redstar_gen_graph($redxml,$outpath."/".basename($redxml).".gen.out"," ");
#   
#   # make dem gen props
#   my $gen_prop = $param->write_harom_hadron_node_ini_xml($props);
#   &run_harom_ptx($gen_prop, $outpath."/".basename($gen_prop).".out"," ");  
#   
#   
#   
#   # copy back 
#   my $file = $param->generalized_perambulator_file(); 
#   my $location = $param->cache_dir() . "/" . $param->stem() . "/gen_props/gen_prop_dbs/dt${dt}/";
#   my $filename = $location . $file;
#   my $scratch_file = $param->scratch_seq_callback("unsmeared_hadron_node_sdb"); 
#   my $cmd = "rcp $scratch_file $filename"; 
#   print "Attempting to execute command: $cmd \n"; 
#   system ( " $cmd " ) == 0 || die ( "failure to rcp : $filename " ) ; 
#   
#   # I hope to someday see this..
#   print "YOU WIN!\n";


# 
# Hopefully this is the new hottness
#


#   # need to do some redstar_gen_graph here to get the smeared/unsmeared list
my $redxml = $param->write_redstar_xml(); 
my $haromxml = $param->write_harom_hadron_node_ini_xml($param->prop_dbs()); 
my $mesonxml = $param->write_meson_hadron_node_ini_xml();

&run_redstar_gen_graph($redxml,$outpath."/".basename($redxml).".gen.out"," ");
&run_harom_cpu($haromxml,$outpath."/".basename($haromxml).".gen_prop.harom.out"," "); 
&run_hadron_node_colorvec($mesonxml,$outpath."/".basename($mesonxml).".hadron_node.out",""); 
&run_redstar_npt($redxml,$outpath."/".basename($redxml).".npt.out"," ");

my $scratch_node = $param->scratch_seq_callback("unsmeared_hadron_node_sdb"); 
my $dest_node = $param->gen_prop_db(); 
my $scratch_corr = $param->scratch_seq_callback("output_sdb"); 
my $dest_corr = $param->work_dir() . "/" . $param->stem(); 
$dest_corr .= "/meson_3pt_redstar/unsmeared_insertion/";
$dest_corr .= "pion_proj0__pion_proj0/";
$dest_corr .= $param->stem() . ".sdb" . $param->seqno(); 

$param->copy_back_rename_rcp($scratch_node,$dest_node);
$param->copy_back_rename_rcp($scratch_corr,$dest_corr);

exit ( 0 ) ; 




#
#
#   RUN CODE
#
#



sub ls_scratch 
{
  my $foo =   `/bin/ls -alh /scratch`;
  print " THIS IS SCRATCH \n $foo \n"; 
}


sub omp_info
{
  my $num_thread = `grep '^processor' /proc/cpuinfo | wc -l ` ;
  chomp $num_thread; 
  $num_thread -= 2; 

  my $nodefile = $ENV{'PBS_NODEFILE'};

  my $num_nodes = `cat $nodefile | wc -l` ;
  chomp $num_nodes; 

  my $run = "env";
  $run .= " OMP_NUM_THREADS=$num_thread";
  $run .= " OMP_PROC_BIND=true"; 
  my $gomp = "0";
  for my $c (1 .. ${num_thread} -1)
  {
    $gomp .= " $c"; 
  }
  $run .= " GOMP_CPU_AFFINITY=\'${gomp}\'";

  print " MPI RUN INFO \n";
  print " num thread = $num_thread \n ";
  print " run = $run \n"; 

  my @r = (); 
  push @r , $run , $num_thread; 
  return \@r ; 
}


sub mpi_info
{
  my $num_thread = `grep '^processor' /proc/cpuinfo | wc -l ` ;
  chomp $num_thread; 
  $num_thread /= 2; 

  my $nodefile = $ENV{'PBS_NODEFILE'};

  my $num_nodes = `cat $nodefile | wc -l` ;
  chomp $num_nodes; 

#  my $base = dirname($0);
#  my $numa = "${base}/numa_script.sh";
  my $numa = $ENV{'HOME'};
  chomp $numa; 
  $numa .= "/scripts/numa_script.sh"; 


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



  print " \n\n FINISHED CHROMA WITH \n\n $chroma_command \n\n";
  &ls_scratch(); 
}


sub run_harom_ptx
{
  my ($input_file , $output_file, $chroma_opts) = @_; 
  my $chroma_exe = "/home/shultz/git-builds/harom-jitptx/gpu-parscalar-Nd3/bin/harom";

  my ($run,$num_nodes) = @{ &mpi_info() }; 

  $chroma_opts .= " -iogeom 1 1 $num_nodes -geom 1 1 $num_nodes "; 

  my $chroma_command = "$run $chroma_exe $chroma_opts -i $input_file -o $output_file ";
 
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($chroma_exe,"run_harom_ptx");

  print " \n\n RUNNING HAROM WITH \n\n $chroma_command \n\n";

  system ( " $chroma_command " ) == 0 || die (" some error for harom, output : $output_file \n $_ ");
  print " \n\n FINISHED HAROM WITH \n\n $chroma_command \n\n";
  &ls_scratch(); 
}

sub run_harom_cpu
{
  my ($input_file , $output_file, $chroma_opts) = @_; 
  my $chroma_exe = "/home/shultz/git-builds/harom/parscalar-Nd3/bin/harom";

  my ($run,$num_nodes) = @{ &omp_info() }; 

  my $chroma_command = "$run $chroma_exe $chroma_opts -i $input_file -o $output_file ";
 
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($chroma_exe,"run_harom_ptx");

  print " \n\n RUNNING HAROM WITH \n\n $chroma_command \n\n";

  system ( " $chroma_command " ) == 0 || die (" some error for harom, output : $output_file \n $_ ");
  print " \n\n FINISHED HAROM WITH \n\n $chroma_command \n\n";
  &ls_scratch(); 
}

sub run_redstar_npt
{
  my ($input_file , $output_file) = @_;
  my $exe = "/home/shultz/git-builds/redstar/bin/redstar_npt";
  my ($run,$num_nodes) = @{ &omp_info() }; 
  my $cmd = "$run $exe $input_file $output_file";
  
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($exe,"run_redstar_npt");

  system ( " $cmd " ) == 0 || die (" some error for redstar_npt, output : $output_file " ); 

  print "FINISHED REDSTAR_NPT \n\n "; 
  &ls_scratch(); 

}

sub run_redstar_gen_graph
{

  my ($input_file , $output_file, $options) = @_;
  my $exe = "/home/shultz/git-builds/redstar/bin/redstar_gen_graph";
  my ($run,$num_nodes) = @{ &omp_info() }; 
  my $cmd = "$run $exe $input_file $output_file ";
  
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }

  &check_exe($exe,"run_redstar_gen_graph");

  system ( " $cmd " ) == 0 || die (" some error for redstar_npt, output : $output_file " ); 

  print "FINISHED REDSTAR GEN GRAPH \n\n";
  &ls_scratch();
}

sub run_hadron_node_colorvec
{ 
  my ($input,$output,$options) = @_; 
  my $exe = "/home/shultz/git-builds/colorvec/bin/hadron_node";
  my $cmd = "$exe $input $output";
  if( -f $output)
  {
    unlink $output; 
  }

  &check_exe($exe,"hadron_node"); 
  
  system(" $cmd " )  == 0 || die("some error for hadron_node, output: $output "); 
  print "FINISHED COLORVEC \n\n";
  &ls_scratch(); 

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








