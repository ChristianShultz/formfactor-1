#!/usr/bin/perl
#
# this file is suppose to be smart enough so I don't have to work hard
#  
#  at the time of writing it runs harom/chroma to get solution vectors/perambulators
#

use strict; 
use Params; 

die "usage: $0 <seqno> <output_path> " unless $#ARGV == 1; 

my $seqno = $ARGV[0];
my $outpath = $ARGV[1];

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
$param->delta_t([28,32,36]);

# configuration stuff
$param->config_type("SCIDAC");
$param->stem("szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per");

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

&run_harom_ptx($src, $outpath."/".$src , " "); 
&run_chroma_ptx($prop, $outpath."/".$prop , " "); 
&run_harom_ptx($peram, $outpath."/".$peram , " "); 



#
#  IF PTX ACTUALLY WORKED WE SHOULD BE DOING GEN PROPS HERE,
#     the perl module that this dude is running knows how to 
#     generate redstar xmlinis and harom_hadron_node inis
#
#  -- the hack is to have the shell script that calls this guy 
#     fire off another qsub for an infiniband harom run to tie
#     all this junk together into a gen prop 
#
#       this is stupid and annoying
#



# shove all the copy back junk into an array
my @rcps = (); 
push @rcps , @{$srcs} , @{$props} , @{$perams};

# copy lots of junk back to cache
print "COPY BACK FILES\n";
foreach my $ele (@rcps)
{
  print $ele . "\n";
}

$param->copy_back_scratch_rcp_arr(\@rcps);

# clean up ?

# code -i inp -o out
#
#    delete inp and out 
#
unlink ( $outpath."/".$src ); 
unlink ( $outpath."/".$prop );
unlink ( $outpath."/".$peram );

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




sub run_chroma_ptx
{
  my ($input_file , $output_file, $chroma_opts) = @_;
  my $chroma_exe = "/home/shultz/git-builds/chroma-jitptx/parscalar-Nd4/bin/chroma";
  my $chroma_command = "$chroma_exe $chroma_opts -i $input_file -o $output_file";
 
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }
 
  system ( " $chroma_command " ) == 0 || die (" some error for chroma, output : $output_file ");
}


sub run_harom_ptx
{
  my ($input_file , $output_file, $chroma_opts) = @_; 
  my $chroma_exe = "/home/shultz/git-builds/harom-jitptx/gpu-parscalar-Nd3/bin/harom";
  my $chroma_command = "$chroma_exe $chroma_opts -i $input_file -o $output_file";
 
  if ( -f $output_file ) 
  {
    unlink $output_file; 
  }
 
  system ( " $chroma_command " ) == 0 || die (" some error for harom, output : $output_file ");
}













