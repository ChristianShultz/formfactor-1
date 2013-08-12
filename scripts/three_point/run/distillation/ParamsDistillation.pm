#!/usr/bin/perl

# make a params package since perl w/o using strict is scary 
package ParamsDistillation;
use File::Copy;
use File::Basename; 
use strict;


# a cute little thing I found on the internet

sub __func__ {(caller(1))[3].'()'}

sub err_func
{
  my ($f,$l,$e) = @_; 
  print "ERROR: $f, line: $l\n  message: $e\n";
  die(); 
}

sub test_func
{
  &err_func(__func__ , __LINE__ , "string");
}


sub new
{
  my $invocant = shift;
  my $class = ref($invocant) || $invocant; 
  my $self =  {
    STEM => "szscl21_20_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per",
    BASEDIR => undef, 
    REDSTARVERSION => 2,
    HAROMVERSION => 1,
    DIAGNOSTIC => 10,
    ZEROUNSMEAREDGRAPHS => "true",
    SEQNO => undef, 
    NUMVECS => 128,
    NTCORR => 40,
    TSOURCE => [0],
    DT => [6],
    PROJOPXML => [],
    L_T =>128,
    L_S => 20,
    DECAYDIR => 3,
    CONVERTUDTOL => undef,
    CONVERTUDTOS => "true", 
    UMASS => "0.05",
    DMASS => "0.05",
    SMASS => "0.1",    
    LINKSMEARTYPE => "STOUT_SMEAR", 
    LINKSMEARFAC => 0.1,
    LINKSMEARNUM => 10,  
    LUSTREDIR => "/lustre/volatile/Spectrum/Clover/NF2+1",
    CACHEDIR => "/lustre/cache/Spectrum/Clover/NF2+1",
    WORKDIR => "/work/JLabLQCD/LHPC/Spectrum/Clover/NF2+1",
    LHPCDIR => "/lustre/cache/LHPC/Spectrum/Clover/NF2+1",
    GRAPHSDB => "graph.sdb",
    OUTPUTSDB => "output.npt.sdb", 
    SMEAREDHADRONNODESDB => "smeared_hadron_node.sdb",
    UNSMEAREDHADRONNODESDB => "unsmeared_hadron_node.sdb",
    SMEAREDHADRONNODEXML => "smeared_hadron_node.xml",
    UNSMEAREDHADRONNODEXML => "unsmeared_hadron_node.xml",
    NUKELIST => "graph.nuke_list.xml",
    REDSTARXML=> "redstar.ini.xml",
    MESONHADRONNODE => "meson_hadron_node.ini.xml",
    HAROMHADRONNODE => "harom_hadron_node.ini.xml",
    DISTSRCINI => "dist_src.ini.xml",
    PROPDISTINI => "prop_dist.ini.xml",
    NPOINTLIST => "npt.list.xml",
    CLOVERCS => 1.589327,
    CLOVERCT => 0.902784,
    XI => 3.5,
    XI_0 => 4.3,
    NU => 1.265,
    RHO => 0.14,
    N_SMEAR => 2, 
    CFGTYPE => "SCIDAC",
    INVERTER => "BICGSTAB", 
    SCRATCHY=>undef,
    SCRATCHYPROPS=>[],
    SCRATCHYPERAMS=>[],
    VERBOSE=>undef, 
    @_, #Override previous attrubutes 
  }; ## self

  return bless $self , $class; 
}


##############################################
## methods to access per-object data        ##
##                                          ##
## With args, they set the value.  Without  ##
## any, they only retrieve it/them.         ##
##############################################

sub stem {
  my $self = shift; 
  if(@_) {$self->{STEM} = shift}
  return $self->{STEM};
}

sub base_dir{
  my $self = shift; 
  if(@_) {$self->{BASEDIR} = shift}
  return $self->{BASEDIR}
}

sub redstar_version{
  my $self = shift; 
  if(@_) {$self->{REDSTARVERSION} = shift}
  return $self->{REDSTARVERSION};
}

sub harom_version{
  my $self = shift; 
  if(@_) {$self->{HAROMVERSION} = shift}
  return $self->{HAROMVERSION};
}

sub diagnostic_level{
  my $self = shift;
  if(@_) {$self->{DIAGNOSTIC} = shift}
  return $self->{DIAGNOSTIC};
}

sub seqno {
  my $self = shift; 
  if(@_) {$self->{SEQNO} = shift}
  return $self->{SEQNO};
}

sub num_vecs {
  my $self = shift; 
  if(@_) {$self->{NUMVECS} = shift}
  return $self->{NUMVECS};
}

sub nt_corr {
  my $self = shift; 
  if(@_) {$self->{NTCORR} = shift}
  return $self->{NTCORR};
}

sub t_sources {
  my $self = shift; 
  if(@_) {$self->{TSOURCE} = shift}
  return $self->{TSOURCE};
}

sub delta_t {
  my $self = shift; 
  if(@_) {$self->{DT} = shift}
  return $self->{DT};
}

sub proj_op_xml {
  my $self = shift; 
  if(@_) {$self->{PROJOPXML} = shift}
  return $self->{PROJOPXML};
}

sub L_t {
  my $self = shift; 
  if(@_) {$self->{L_T} = shift}
  return $self->{L_T}
}

sub L_s {
  my $self = shift; 
  if(@_) {$self->{L_S} = shift}
  return $self->{L_S}
}

sub decay_dir {
  my $self = shift; 
  if(@_) {$self->{DECAYDIR} = shift}
  return $self->{DECAYDIR}
}

sub convertUDtoL {
  my $self = shift; 
  if(@_) {$self->{CONVERTUDTOL} = shift}
  return $self->{CONVERTUDTOL}
}

sub convertUDtoS {
  my $self = shift; 
  if (@_) {$self->{CONVERTUDTOS} = shift}
  return $self->{CONVERTUDTOS}
}


sub u_mass {
  my $self = shift; 
  if(@_) {$self->{UMASS} = shift}
  return $self->{UMASS};
}

sub d_mass {
  my $self = shift; 
  if(@_) {$self->{DMASS} = shift}
  return $self->{DMASS};
}

sub s_mass {
  my $self = shift; 
  if(@_) {$self->{SMASS} = shift}
  return $self->{SMASS};
}

sub link_smear_type{
  my $self = shift; 
  if(@_) {$self->{LINKSMEARTYPE} = shift}
  return $self->{LINKSMEARTYPE};
}

sub link_smear_fac{
  my $self = shift; 
  if(@_) {$self->{LINKSMEARFAC} = shift}
  return $self->{LINKSMEARFAC};
}

sub link_smear_num{
  my $self = shift; 
  if(@_) {$self->{LINKSMEARNUM} = shift}
  return $self->{LINKSMEARNUM};
}

sub lustre_dir {
  my $self = shift; 
  if(@_) {$self->{LUSTREDIR} = shift}
  return $self->{LUSTREDIR};
}

sub cache_dir {
  my $self = shift; 
  if(@_) {$self->{CACHEDIR} = shift}
  return $self->{CACHEDIR};
}

sub work_dir {
  my $self = shift; 
  if(@_) {$self->{WORKDIR} = shift}
  return $self->{WORKDIR};
}

sub lhpc_dir {
  my $self = shift;
  if(@_) {$self->{LHPCDIR} = shift}
  return $self->{LHPCDIR};
}

sub graph_sdb {
  my $self = shift; 
  if(@_) {$self->{GRAPHSDB} = shift}
  return $self->{GRAPHSDB};
}

sub output_sdb {
  my $self = shift; 
  if(@_) {$self->{OUTPUTSDB} = shift}
  return $self->{OUTPUTSDB};
}

sub smeared_hadron_node_sdb
{
  my $self = shift;
  if(@_) {$self->{SMEAREDHADRONNODESDB} = shift}
  return $self->{SMEAREDHADRONNODESDB};
}

sub unsmeared_hadron_node_sdb
{
  my $self = shift;
  if(@_) {$self->{UNSMEAREDHADRONNODESDB} = shift}
  return $self->{UNSMEAREDHADRONNODESDB};
}

sub smeared_hadron_node_xml
{
  my $self = shift;
  if(@_) {$self->{SMEAREDHADRONNODEXML} = shift}
  return $self->{SMEAREDHADRONNODEXML};
}

sub unsmeared_hadron_node_xml
{
  my $self = shift;
  if(@_) {$self->{UNSMEAREDHADRONNODEXML} = shift}
  return $self->{UNSMEAREDHADRONNODEXML};
}

sub graph_nuke_list_xml
{
  my $self = shift;
  if(@_) {$self->{NUKELIST} = shift}
  return $self->{NUKELIST};
}

sub redstar_ini {
  my $self = shift; 
  if(@_) {$self->{REDSTARXML} = shift}
  return $self->{REDSTARXML};
}

sub meson_hadron_node {
  my $self = shift; 
  if(@_) {$self->{MESONHADRONNODE} = shift}
  return $self->{MESONHADRONNODE};
}

sub harom_hadron_node {
  my $self = shift; 
  if(@_) {$self->{HAROMHADRONNODE} = shift}
  return $self->{HAROMHADRONNODE};
}

sub dist_src_ini {
  my $self = shift; 
  if(@_) {$self->{DISTSRCINI} = shift}
  return $self->{DISTSRCINI};
}

sub prop_dist_ini{
  my $self = shift;
  if(@_) {$self->{PROPDISTINI} = shift}
  return $self->{PROPDISTINI};
}

sub peram_dist_ini{
  my $self = shift;
  if(@_) {$self->{PERAMDISTINI} = shift}
  return $self->{PERAMDISTINI};
}

sub npt_list_xml{
  my $self = shift;
  if(@_) {$self->{NPOINTLIST} = shift}
  return $self->{NPOINTLIST};
}

sub clover_c_s{
  my $self = shift; 
  if(@_) {$self->{CLOVERCS} = shift}
  return $self->{CLOVERCS}
}

sub clover_c_t{
  my $self = shift; 
  if(@_) {$self->{CLOVERCT} = shift}
  return $self->{CLOVERCT}
}

sub xi{
  my $self = shift; 
  if(@_) {$self->{XI} = shift}
  return $self->{XI}
}

sub xi_0{
  my $self = shift; 
  if(@_) {$self->{XI_0} = shift}
  return $self->{XI_0}
}

sub nu{
  my $self = shift; 
  if(@_) {$self->{NU} = shift}
  return $self->{NU}
}

sub rho{
  my $self = shift; 
  if(@_) {$self->{RHO} = shift}
  return $self->{RHO}
}

sub n_smear{
  my $self = shift; 
  if(@_) {$self->{N_SMEAR} = shift}
  return $self->{N_SMEAR}
}

sub config_type{
  my $self = shift; 
  if(@_) {$self->{CFGTYPE} = shift}
  return $self->{CFGTYPE}
}

sub inverter{
  my $self = shift; 
  if(@_) {$self->{INVERTER} = shift}
  return $self->{INVERTER};
}

sub scratchy{
  my $self = shift; 
  if(@_) {$self->{SCRATCHY} = shift}
  return $self->{SCRATCHY};
}

sub scratchy_props{
  my $self = shift; 
  if (@_) {$self->{SCRATCHYPROPS} = shift}
  return $self->{SCRATCHYPROPS};
}

sub scratchy_perams{
  my $self = shift; 
  if (@_) {$self->{SCRATCHYPERAMS} = shift}
  return $self->{SCRATCHYPERAMS};
}

sub verbose{
  my $self = shift; 
  if (@_) {$self->{VERBOSE} = shift}
  return $self->{VERBOSE};
}

#
# More methods.. no args
#

# Robert's randomized origin..
sub t_origin{

  my $self = shift;

  srand($self->seqno());

  foreach my $i (1..20)
  {
    rand (1.0);
  }

  return int(rand($self->L_t())); 
}

# to avoid calling stem 
sub ensemble{
  my $self = shift; 
  return $self->{STEM};
}


# one sub to rule them all

sub mass_string {

  my $self = shift; 
  my $mass = shift; 

  if ($mass == $self->u_mass())
  {
    $mass = "U" . $mass;
  }
  elsif($mass == $self->s_mass())
  {
    $mass = "S" . $mass; 
  }
  else
  {
    die("unrecognized quark mass $mass \n");
  }

  return $mass;
}

# copy a file to scratch and return the string of its location

sub copy_file_rename
{
  my $self = shift; 
  my $full = shift;
  my $dest = shift;  

  print "In copy_file_to_scratch, copying $full to $dest \n" unless ! $self->verbose();
  if(! -e $dest)
  {
    copy($full, $dest) || die("Couldn't copy $full\n to $dest\n $_");
  }
  else
  {
    print "skipped copying $dest since a copy already existed\n";
  }
}

sub copy_file_to_scratch
{
  my $self = shift; 
  my $full = shift; 
  my $file = basename($full);

  print "In copy_file_to_scratch, copying $file \n" unless ! $self->verbose();
  if(! -e "/scratch/${file}")
  {
    copy($full,"/scratch/${file}") || die("Couldn't copy $file \nfull: $full \n\n $_");
  }
  else
  {
    print "skipped copying $file since a copy already existed at /scratch/${file}\n";
  }
  return "/scratch/${file}";
}

# scratchify a file name

sub scratchify_file_name
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1; 
  my $self = shift;
  my $full = shift; 
  my $file = basename($full);
  return "/scratch/$file";
}

sub scratchify_file_name_arr
{
  my $self = shift; 
  my $ref = shift; 
  my @a = @{ $ref } ;
  my @r = (); 
  foreach my $ele (@a)
  {
    push @r, $self->scratchify_file_name($ele); 
  }

  return \@r;
}

sub copy_back_rename_rcp
{
  my $self = shift; 
  my $file = shift; 
  my $dest = shift; 
  my $cmd = "rcp $file $dest";
  print __func__ . " trying to execute command\n $cmd\n "; 
  system ( " $cmd " ) == 0 || die ( " failure to rcp : $dest \n $_ " ); 
}

sub copy_back_scratch_rcp
{
  my $self = shift; 
  my $dest = shift; 
  my $file = basename($dest); 
  my $cmd = "rcp /scratch/${file} $dest";
  print "copy_back_scratch_rcp : trying to execute command\n $cmd\n "; 
  system ( " $cmd " ) == 0 || die ( " failure to rcp : $dest \n $_ " ); 
}

sub copy_back_scratch_rcp_arr
{
  my $self = shift; 
  my $ref = shift; 
  my @a = @{ $ref };

  foreach my $ele (@a) 
  {
    $self->copy_back_scratch_rcp($ele); 
  }
}



#
#
#  SDBS THAT ARE ALREADY ON DISK
#
#

sub gen_prop_db
{ 
  my $self = shift; 
  my $f = $self->cache_dir() . "/" . $self->stem();
  $f .= "/gen_props/gen_prop_dbs/dt"; 
  my @dts = @{$self->delta_t()}; 
  die ("confused") unless $#dts == 0; 
  $f .= $dts[0] . "/" . $self->stem() . ".unsmeared_nodes.sdb" . $self->seqno();
  return $f;
}


sub prop_dbs{
  my $self = shift; 

#    if ( $self->scratchy() )
#    {
#  
#      if ( $self->verbose() ) 
#      {
#        my @t = @{$self->scratchy_props()};
#        print "prop_dbs props @t \n";
#      }
#  
#      return $self->scratchy_props(); 
#    } 


  my @prop_dbs; 
#    my $seqno = $self->seqno();
#    my $t_ref = $self->t_sources(); 
#    my @flav = ("light");
#    foreach my $tt (@$t_ref)
#    {
#      foreach my $flavor (@flav)
#      {
#        push @prop_dbs , "$cache_dir/${stem}/prop_db/${stem}.prop.n128.${flavor}.t0_${tt}.sdb${seqno}";
#      }
#    }
#  
# don't need diag guys since no bubbles


  my $loc = $self->cache_dir() . "/" . $self->stem() . "/prop_mod/";
  my $f = $loc . $self->stem() . ".prop.n" . $self->num_vecs(); 
  $f .= ".t0_";
  my $tail = ".mod" . $self->seqno(); 

  push @prop_dbs , $f . "0" . $tail ;
  
  foreach my $t (@{$self->delta_t()})
  {
    push @prop_dbs , $f . $t . $tail ; 
  }


  my @scratch_dbs = ();

  foreach my $db (@prop_dbs)
  {
    push @scratch_dbs , $self->copy_file_to_scratch($db);
  }

# push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_0-124_inc4.sdb${seqno}");
# push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_1-125_inc4.sdb${seqno}");
# push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_2-126_inc4.sdb${seqno}");
# push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_3-127_inc4.sdb${seqno}");

  return \@scratch_dbs; 
}

sub peram_dbs
{
  my $self = shift; 

  my $loc = $self->cache_dir() . "/" . $self->stem() . "/prop_db/";
  my $f = $loc . $self->stem() . ".prop.n" . $self->num_vecs(); 
  $f .= ".t0_";
  my $tail = ".sdb" . $self->seqno(); 
  my @prop_dbs = (); 
  push @prop_dbs , $f . "0" . $tail ;
  
  foreach my $t (@{$self->delta_t()})
  {
    push @prop_dbs , $f . $t . $tail ; 
  }

  my @scratch_dbs = ();
  foreach my $db (@prop_dbs)
  {
    push @scratch_dbs , $self->copy_file_to_scratch($db);
  }

  return \@scratch_dbs; 
}


sub meson_dbs
{
  my $self = shift;
  my @mdbs = ();

  my $stem = $self->stem(); 
  my $cache_dir = $self->cache_dir(); 
  my $cache_stem = "$cache_dir/${stem}/meson_db_dispmom/${stem}.meson.colorvec.";
  push @mdbs , $cache_stem . "1disp.sdb" . $self->seqno();
  push @mdbs , $cache_stem . "2disp.sdb" . $self->seqno();
 
  # there seem to be multiple copies of this guy?? -- use only the cache guys
  # my $lhpc_dir = $self->lhpc_dir(); 
  # my $lhpc_stem = "$lhpc_dir/${stem}/meson_db_combined/${stem}.meson.colorvec.sdb";
  # push @mdbs , $lhpc_stem . $self->seqno(); 

  my $cache_stem2 = "$cache_dir/${stem}/meson_db_combined/${stem}.meson.colorvec.sdb";
  push @mdbs , $cache_stem2 . $self->seqno(); 

  my @scratch_dbs = ();

  foreach my $db (@mdbs)
  {
    push @scratch_dbs , $self->copy_file_to_scratch($db); 
  }
  return \@scratch_dbs; 
}

sub baryon_dbs
{
  my $self = shift; 
  my @bdbs;
  return \@bdbs; 
}


sub glue_dbs
{
  my $self = shift; 
  my @gdbs;
  return \@gdbs; 
}

sub tetraquark_dbs
{
  my $self = shift; 
  my @tet;
  return \@tet; 
}


sub sequenceify 
{
  my $self = shift; 
  my $name = shift; 

  $name = $name . $self->seqno(); 

  return $name; 
}

sub scratch_and_sequenceify
{
  my $self = shift; 
  my $x = shift;  
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}


# basically do some callb(l)ack magic

sub scratch_seq_callback
{
  my ( $self, $string ) = @_;

  my %hash = (); 

  # SDB NAMES
  $hash{"smeared_hadron_node_sdb"} = \&ParamsDistillation::smeared_hadron_node_sdb;
  $hash{"unsmeared_hadron_node_sdb"} = \&ParamsDistillation::unsmeared_hadron_node_sdb;
  $hash{"graph_sdb"} = \&ParamsDistillation::graph_sdb;
  $hash{"output_sdb"} = \&ParamsDistillation::output_sdb;

  # XML NAMES
  $hash{"smeared_hadron_node_xml"} = \&ParamsDistillation::smeared_hadron_node_xml;
  $hash{"unsmeared_hadron_node_xml"} = \&ParamsDistillation::unsmeared_hadron_node_xml;

  # REDSTAR
  $hash{"redstar_ini"} = \&ParamsDistillation::redstar_ini;
  $hash{"nuke_list"} = \&ParamsDistillation::graph_nuke_list_xml; 

  # COLORVEC
  $hash{"meson_hadron_node"} = \&ParamsDistillation::meson_hadron_node; 

  # HAROM
  $hash{"harom_hadron_node"} = \&ParamsDistillation::harom_hadron_node;
  $hash{"peram_dist_ini"} = \&ParamsDistillation::peram_dist_ini;
  $hash{"dist_src_ini"} = \&ParamsDistillation::dist_src_ini;

  # CHROMA
  $hash{"prop_dist_ini"} = \&ParamsDistillation::prop_dist_ini;

  if ( exists $hash{$string} ) 
  {
    # black magic! -- store a ref to a sub in the hash then call it 
    #     using the "self" as input so we keep all of our info then 
    #     pass the result (some filename string)it to the 
    #     scratchification/sequenceification function

    print __func__ . ": running $string \n" unless ! $self->verbose() ;
    return $self->scratch_and_sequenceify( &{ $hash{$string} }($self) ); 
  }
  else
  {
    my @k = keys %hash; 
    print "ERROR: don't know what a $string is, a smarter set of choices is\n";
    foreach my $kk (@k) 
    {
      print $kk . "\n";
    }
    die ( "unholy death" );
  }
}

sub stupid
{
  my $self = shift; 
  print "stupid \n";
}

#
#
#
#
#   HAROM AND CHROMA STUFF FOR GEN PROP GENERATION
#
#
#
#
#



#
#  FILENAME/LOCATION UTILITY FUNCTIONS
#

sub distillation_inversion_source_file
{
  my ($self,$f_or_b) = @_; 

  my $f = undef; 
  $f = "srcf.mod" unless $f_or_b ne "f"; 
  $f = "srcb.mod" unless $f_or_b ne "b"; 
  &err_func(__func__,__LINE__,"wanted f or b, got $f_or_b") unless $f; 
  return $self->scratch_and_sequenceify($f); 
}

sub distillation_propagator_file
{
  my ($self,$f_or_b) = @_; 

  my $f = undef; 
  $f = "propf.mod" unless $f_or_b ne "f"; 
  $f = "propb.mod" unless $f_or_b ne "b"; 
  &err_func(__func__,__LINE__,"wanted f or b, got $f_or_b") unless $f; 
  return $self->scratch_and_sequenceify($f); 
}

sub generalized_perambulator_file
{
  my $self = shift; 
  my $dt = $self->arr_string("_",$self->delta_t());
  my $ts = $self->arr_string("_",$self->t_sources()); 
  my $nv = $self->num_vecs(); 
  my $f = $self->stem() . ".n" . $nv . ".t0" . $ts . ".dt" . $dt . ".sdb"; 
  return $self->sequenceify($f); 
}

sub colorvec_file
{
  my $self = shift; 
  my $c = ""; 
  my $cf = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->stem() . "/eigs_mod/" . $self->stem() . ".3d.eigs.mod" . $self->seqno(); 
  $c = $self->copy_file_to_scratch($cf); 
  return $c; 
}


sub gauge_file_lime
{
  my $self = shift; 
  my $g = ""; 
  my $path = "/lustre/cache/LHPC/Spectrum/Clover/NF2+1/" . $self->stem() . "/cfgs/";
  my $gf = $path . $self->stem() . "_cfg_" . $self->seqno() . ".lime"; 
  $g = $self->copy_file_to_scratch($gf); 

  return $g;
}

sub gauge_file_3d
{
  my $self = shift; 
  my $g = "";

  my $gfpath = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->stem() . "/cfgs_mod/";
  my $gf = $gfpath . $self->stem() . ".3d.gauge.mod" . $self->seqno(); 
  $g = $self->copy_file_to_scratch($gf); 
  return $g; 
}


sub print_elem_arr
{
  my $self = shift; 
  my $ref = shift; 
  my @a = @{ $ref };

  foreach my $elem (@a)
  {
    print OUT <<EOF;
      <elem>$elem</elem>
EOF
  }
}

sub arr_string
{
  my $self = shift; 
  my $delim = shift; 
  my $ref = shift; 
  my @a = @{ $ref };  
  my $r = "";

  foreach my $ele (@a)
  {
    $r .= $delim . $ele; 
  }

  return $r; 
}

sub time_source_string
{
  &err_func(__func__,__LINE__,"expected 2 args") unless $#_ == 2; 
  my $self = shift;
  my $delim = shift; 
  my $line_type = shift;  


  print __func__ . " delim = 1 $delim 1 , line = $line_type \n" unless ! $self->verbose();

  my @t = (); 

  if($line_type eq "f")
  {
    @t = @{$self->t_sources()};
  }
  elsif ($line_type eq "b")
  {
    my @tmp = @{ $self->t_sources() };
    my @tmp2 = @{ $self->delta_t() };

    foreach my $t1 (@tmp)
    {
      foreach my $t2 (@tmp2)
      {
        push @t , $t1 + $t2; 
      }
    }
  }
  elsif ($line_type eq "a")
  {
    # do nothing
  }
  else
  {
    print "ERROR: unknown quark line type, $line_type, exiting.. \n";
    die("unholy death");  
  }

  return $self->arr_string($delim,\@t); 
}



#
#
#
#
#
#  XML UTILITY  
#
#
#
#
#



sub print_harom_header
{
  my $self = shift; 

  print OUT <<EOF;
<?xml version="1.0"?>
<harom>
<Param>
  <InlineMeasurements>
EOF

}

sub print_harom_trailer
{
  my $self = shift; 

  print OUT <<EOF;
  </InlineMeasurements>
  <nrow>$self->{L_S} $self->{L_S} $self->{L_S} $self->{L_T}</nrow>
</Param>

<RNG>
  <Seed>
    <elem>11</elem>
    <elem>11</elem>
    <elem>11</elem>
    <elem>0</elem>
  </Seed>
</RNG>

</harom>
EOF

}


sub print_chroma_header
{
  my $self = shift; 

  print OUT <<EOF;
<?xml version="1.0"?>
<chroma>
<Param>
  <InlineMeasurements>
EOF

}

sub print_chroma_trailer
{
  my $self = shift; 
  my $gauge_file = $self->gauge_file_lime(); 

  print OUT <<EOF;
  </InlineMeasurements>
  <nrow>$self->{L_S} $self->{L_S} $self->{L_S} $self->{L_T}</nrow>
</Param>

<RNG>
  <Seed>
    <elem>11</elem>
    <elem>11</elem>
    <elem>11</elem>
    <elem>0</elem>
  </Seed>
</RNG>

<Cfg>
  <cfg_type>$self->{CFGTYPE}</cfg_type>
  <cfg_file>$gauge_file</cfg_file>
</Cfg>

</chroma>
EOF

}

#
#
#  DISTILLATION 
#
#


# harom solution vector sources in distillation

sub print_MAKE_SOURCE_DISTILLATION__harom
{
  &err_func(__func__,__LINE__,"expected 3 args") unless $#_ == 3; 

  my ($self,$tsource,$colorvec,$src) = @_; 
  my $nv = $self->num_vecs();  

  print OUT<<EOF;

    <elem>
      <annotation>
        construct propagator source vectors for chroma
      </annotation>
      <Name>MAKE_SOURCE_DISTILLATION</Name>
      <Frequency>1</Frequency>
      <Param>
        <Contractions>
          <num_vecs>$nv</num_vecs>
          <t_sources>$tsource</t_sources>
        </Contractions>
      </Param>
      <NamedObject>
        <colorvec_file>$colorvec</colorvec_file>
        <src_file>$src</src_file>
      </NamedObject>
    </elem>

EOF

}


sub print_PROP_DISTILLATION_contraction_xml__chroma
{
  &err_func(__func__,__LINE__,"expected 4 args") unless $#_ == 4;
  my ($self,$mass,$tsource,$ntf,$ntb) = @_; 
  my $nv = $self->num_vecs();

  my $m = $self->mass_string($mass); 

  print OUT<<EOF;
        <Contractions>
          <mass_label>$m</mass_label>
          <num_vecs>$nv</num_vecs>
          <t_sources>$tsource</t_sources>
          <decay_dir>3</decay_dir>
          <Nt_forward>$ntf</Nt_forward>
          <Nt_backward>$ntb</Nt_backward>
        </Contractions>
EOF
}

sub print_PROP_DISTILLATION_clover_propagator_xml__chroma
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1;
  my ($self,$mq) = @_;

  print OUT<<EOF;
        <Propagator>
          <version>10</version>
          <quarkSpinType>FULL</quarkSpinType>
          <obsvP>false</obsvP>
          <numRetries>1</numRetries>
          <FermionAction>
            <FermAct>CLOVER</FermAct>
            <Mass>$mq</Mass>
            <clovCoeffR>$self->{CLOVERCS}</clovCoeffR>
            <clovCoeffT>$self->{CLOVERCT}</clovCoeffT>
            <AnisoParam>
              <anisoP>true</anisoP>
              <t_dir>3</t_dir>
              <xi_0>$self->{XI_0}</xi_0>
              <nu>$self->{NU}</nu>
            </AnisoParam>
            <FermState>
              <Name>STOUT_FERM_STATE</Name>
              <rho>$self->{RHO}</rho>
              <n_smear>$self->{N_SMEAR}</n_smear>
              <orthog_dir>3</orthog_dir>
              <FermionBC>
                <FermBC>SIMPLE_FERMBC</FermBC>
                <boundary>1 1 1 -1</boundary>
              </FermionBC>
            </FermState>
          </FermionAction>
EOF
  $self->print_PROP_DISTILLATION_inverter_xml__chroma($mq);
  print OUT<<EOF;
         </Propagator>
EOF
}

sub print_PROP_DISTILLATION_inverter_xml__chroma
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1;
  my ($self,$mq) = @_; 

  my %hash = ();

  $hash{"BICGSTAB"} = \&ParamsDistillation::print_BICGSTAB_inverter_xml__chroma;
  $hash{"QUDA_CLOVER_INVERTER"} = \&ParamsDistillation::print_QUDA_CLOVER_INVERTER_inverter_xml__chroma;
  if ( exists $hash{$self->inverter()} ) 
  {
    &{ $hash{$self->inverter()} }($self,$mq) ; 
  }
  else
  {
    my @k = keys %hash; 
    print "ERROR: don't know what a ". $self->inverter() ." is, a smarter set of choices is\n";
    foreach my $kk (@k) 
    {
      print $kk . "\n";
    }
    die ( "unholy death" );
  }
}

sub print_BICGSTAB_inverter_xml__chroma
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1;
  my ($self,$mq) = @_; 
  print OUT<<EOF;
          <InvertParam>
            <invType>BICGSTAB_INVERTER</invType>
            <RsdBiCGStab>1e-7</RsdBiCGStab>
            <MaxBiCGStab>1000</MaxBiCGStab>
          </InvertParam>
EOF
}


sub print_QUDA_CLOVER_INVERTER_inverter_xml__chroma
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1;
  my ($self,$mq) = @_; 
  print OUT<<EOF;
        <InvertParam>
          <invType>QUDA_CLOVER_INVERTER</invType>
          <CloverParams>
            <Mass>$mq</Mass>
            <clovCoeffR>$self->{CLOVERCS}</clovCoeffR>
            <clovCoeffT>$self->{CLOVERCT}</clovCoeffT>
            <AnisoParam>
              <anisoP>true</anisoP>
              <t_dir>3</t_dir>
              <xi_0>$self->{XI_0}</xi_0>
              <nu>$self->{NU}</nu>
            </AnisoParam>
          </CloverParams>
          <RsdTarget>5.0e-7</RsdTarget>
          <RsdToleranceFactor>10</RsdToleranceFactor>
          <SilentFailP>false</SilentFailP>
          <MaxIter>10000</MaxIter>
          <AntiPeriodicT>true</AntiPeriodicT>
          <SolverType>BICGSTAB</SolverType>
          <Verbose>false</Verbose>
          <AsymmetricLinop>false</AsymmetricLinop>
          <CudaReconstruct>RECONS_12</CudaReconstruct>
          <!-- CudaSloppyPrecision>SINGLE</CudaSloppyPrecision -->
          <!-- Delta>0.001</Delta -->
          <CudaSloppyPrecision>HALF</CudaSloppyPrecision>
          <Delta>0.1</Delta>
          <CudaSloppyReconstruct>RECONS_12</CudaSloppyReconstruct>
          <AxialGaugeFix>false</AxialGaugeFix>
          <AutotuneDslash>true</AutotuneDslash>
        </InvertParam>
EOF
}

sub print_PROP_DISTILLATION_xml__chroma
{
  &err_func(__func__,__LINE__,"expected 5 args") unless $#_ == 5;
  my ($self,$mass,$f_or_b,$g,$src,$soln) = @_;

  my $tsource = $self->time_source_string(" ",$f_or_b); 
  my $ntf = 0; 
  my $ntb = 0; 

  $ntf = $self->nt_corr() unless ! $f_or_b eq "f"; 
  $ntb = $self->nt_corr() unless ! $f_or_b eq "b";
  if ( ($ntf== 0) and ($ntb == 0) )
  {
    &err_func(__func__,__LINE__,"passed $f_or_b, should be f or b"); 
  }

  print OUT<<EOF;
      <elem>
        <Name>PROP_DISTILLATION</Name>
        <Frequency>1</Frequency>
        <Param>
EOF

  $self->print_PROP_DISTILLATION_contraction_xml__chroma($mass,$tsource,$ntf,$ntb);
  $self->print_PROP_DISTILLATION_clover_propagator_xml__chroma($mass); 

  print OUT<<EOF;
        </Param>
        <NamedObject>
          <gauge_id>$g</gauge_id>
          <src_file>$src</src_file>
          <soln_file>$soln</soln_file>
        </NamedObject>
      </elem>
EOF

}



sub print_link_smearing_xml__generic
{
  my $self = shift; 
  my $link_smear_type = $self->link_smear_type(); 
  my $link_smear_fac = $self->link_smear_fac();
  my $link_smear_num = $self->link_smear_num();

  # hardwire no smear temporal 
  print OUT<<EOF;
        <LinkSmearing>
          <LinkSmearingType>$link_smear_type</LinkSmearingType>
          <link_smear_fact>$link_smear_fac</link_smear_fact>
          <link_smear_num>$link_smear_num</link_smear_num>
          <no_smear_dir>3</no_smear_dir>
        </LinkSmearing>

EOF

}


sub print_flavor_to_mass_xml__generic
{
  my $self = shift; 
  my $um = $self->mass_string($self->u_mass()); 
  my $sm = $self->mass_string($self->s_mass()); 

  print OUT<<EOF;
      <FlavorToMass>
        <elem>
          <flavor>l</flavor>
          <mass>$um</mass>
        </elem>
EOF
#        <elem>
#          <flavor>s</flavor>
#          <mass>$sm</mass>
#        </elem>
  print OUT<<EOF;
      </FlavorToMass>
EOF
}






sub print_UNSMEARED_HADRON_NODE_DISTILLATION__harom
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1; 
  my ($self , $solnfiles_ref ) = @_; 
  my $t_origin = $self->t_origin();

  my $gauge = $self->gauge_file_3d(); 
  my $unhadnode = $self->scratch_seq_callback("unsmeared_hadron_node_xml"); 
  my $unhadop = $self->scratch_seq_callback("unsmeared_hadron_node_sdb");


  print OUT<<EOF;

  <elem>
    <Name>UNSMEARED_HADRON_NODE_DISTILLATION</Name>
    <Frequency>1</Frequency>
    <Param>
      <version>$self->{HAROMVERSION}</version>
      <num_vecs>$self->{NUMVECS}</num_vecs>
      <displacement_length>1</displacement_length>
      <t_origin>$t_origin</t_origin>
EOF

  $self->print_flavor_to_mass_xml__generic();
  $self->print_link_smearing_xml__generic();

  print OUT<<EOF;
    </Param>
    <NamedObject>
      <hadron_node_xmls>
        <elem>$unhadnode</elem>
      </hadron_node_xmls>
        <gauge_file>$gauge</gauge_file>
        <soln_files>
EOF

  $self->print_elem_arr($solnfiles_ref);

  print OUT<<EOF;
        </soln_files>
        <hadron_op_file>$unhadop</hadron_op_file>
    </NamedObject>
  </elem>
EOF


}



#
#
# WRITERS 
#
#

sub write_harom_hadron_node_ini_xml
{
  print __func__ . ": entering\n";
  &err_func(__func__,__LINE__,"expected 1 arg (array ref of soln vectors)") unless $#_ == 1; 
  my ($self,$ref) = @_;  
  my $file = $self->scratch_seq_callback("harom_hadron_node");
  open OUT , ">" , $file;
  print __func__ . ": writing harom header\n";
  $self->print_harom_header();
  print __func__ . ": writing UNSMEARED_HAD_NODE_DIST header\n";
  $self->print_UNSMEARED_HADRON_NODE_DISTILLATION__harom($ref);
  print __func__ . ": writing harom trailer\n";
  $self->print_harom_trailer();
  close OUT;
  print __func__ . ": exitingg\n";
  return $file;  
}

sub write_harom_distillation_source_ini_xml
{
  my $self = shift; 
  my $file = $self->scratch_seq_callback("dist_src_ini");
  my $colorvec = $self->colorvec_file(); 
  my @inversion_sources = (); 
  open OUT , ">" , $file; 
  $self->print_harom_header();

  # forward
  my $tsource = $self->time_source_string(" ","f");
  my $src = $self->distillation_inversion_source_file("f");  
  push @inversion_sources , $src;
  $self->print_MAKE_SOURCE_DISTILLATION__harom($tsource,$colorvec,$src);

  # backward
  my $tsource = $self->time_source_string(" ","b");
  my $src = $self->distillation_inversion_source_file("b");  
  push @inversion_sources, $src; 
  $self->print_MAKE_SOURCE_DISTILLATION__harom($tsource,$colorvec,$src);

  $self->print_harom_trailer();
  close OUT; 

  my @r = (); 
  push @r, $file; 
  push @r, \@inversion_sources; 
  return \@r; 
}

sub write_chroma_distillation_prop_ini_xml
{
  my $self = shift; 
  my $file = $self->scratch_seq_callback("prop_dist_ini");
  my @props = (); 
  open OUT , ">" , $file;
  my $id = "default_gauge_field"; 
  $self->print_chroma_header();

  # forward
  my $src = $self->distillation_inversion_source_file("f");
  my $prop = $self->distillation_propagator_file("f");
  push @props , $prop; 
  $self->print_PROP_DISTILLATION_xml__chroma($self->u_mass(),   
    "f", $id, $src, $prop); 

  # backward
  my $src = $self->distillation_inversion_source_file("b");
  my $prop = $self->distillation_propagator_file("b");
  push @props , $prop; 
  $self->print_PROP_DISTILLATION_xml__chroma($self->u_mass(),   
    "b", $id, $src, $prop); 

  $self->print_chroma_trailer(); 
  close OUT; 
  my @r = (); 
  push @r, $file; 
  push @r, \@props;
  return \@r;
}












#
#
#
#
#   REDSTAR STUFF 
#
#
#
#



sub print_redstar_npoint_list
{
  my $self = shift; 
  my $file = $self->npt_list_xml(); 
  die("$file does not exist") unless (-e $file);

  open IN , "<" , $file; 
  while (<IN>)
  {
    print OUT $_; 
  }
  close IN; 
}

sub print_redstar_Param()
{
  my $self = shift; 
  print OUT <<EOF;
  <Param>
    <version>$self->{REDSTARVERSION}</version>
    <diagnostic_level>$self->{DIAGNOSTIC}</diagnostic_level>
    <Nt_corr>$self->{NTCORR}</Nt_corr>
    <average_1pt_diagrams>true</average_1pt_diagrams>
    <zeroUnsmearedGraphsP>$self->{ZEROUNSMEAREDGRAPHS}</zeroUnsmearedGraphsP>
EOF
  if($self->convertUDtoL())
  {
    print OUT <<EOF;
    <convertUDtoL>$self->{CONVERTUDTOL}</convertUDtoL>
EOF
  }
  if($self->convertUDtoS())
  {
    print OUT <<EOF; 
    <convertUDtoS>$self->{CONVERTUDTOS}</convertUDtoS>
EOF
  }
  print OUT <<EOF;
    <Layout>
      <lattSize>$self->{L_S} $self->{L_S} $self->{L_S} $self->{L_T}</lattSize>
      <decayDir>$self->{DECAYDIR}</decayDir>
    </Layout>
EOF

  $self->print_redstar_npoint_list(); 

  print OUT <<EOF;
  </Param> 
EOF

}

sub print_redstar_DBFiles
{
  my $self = shift;
  my $proj_ops = $self->proj_op_xml(); 
  my @projops = @$proj_ops; 
  my $had_graph = $self->scratch_seq_callback("graph_sdb"); 
  my $shad_node_xml = $self->scratch_seq_callback("smeared_hadron_node_xml"); 
  my $shad_node_out_db = $self->scratch_seq_callback("smeared_hadron_node_sdb");
  my $uhad_node_xml = $self->scratch_seq_callback("unsmeared_hadron_node_xml"); 
  my $uhad_node_out_db = $self->scratch_seq_callback("unsmeared_hadron_node_sdb");
  my $out_db = $self->scratch_seq_callback("output_sdb"); 
  my $nuke_list = $self->scratch_seq_callback("nuke_list");

  print OUT <<EOF;
  <DBFiles>
    <proj_op_xmls>
EOF

  foreach (@projops)
  {
    print OUT <<EOF;
      <elem>$_</elem>
EOF
  }    

  print OUT <<EOF; 
    </proj_op_xmls>
    <noneval_graph_xml>$nuke_list</noneval_graph_xml>
    <smeared_hadron_node_xml>$shad_node_xml</smeared_hadron_node_xml>
    <unsmeared_hadron_node_xml>$uhad_node_xml</unsmeared_hadron_node_xml>
    <hadron_npt_graph_db>$had_graph</hadron_npt_graph_db>
    <hadron_node_dbs>
      <elem>$shad_node_out_db</elem>
      <elem>$uhad_node_out_db</elem>
    </hadron_node_dbs>
    <output_db>$out_db</output_db>
  </DBFiles>
EOF
}

sub write_redstar_xml
{
  my $self = shift;

  my $file = $self->scratch_seq_callback("redstar_ini"); 

  if(-e $file)
  {
    unlink($file) || die($_); 
  }
  open OUT , ">" , $file;

  print OUT <<EOF;
<?xml version="1.0"?>
<RedstarNPt>
EOF
  $self->print_redstar_Param();
  $self->print_redstar_DBFiles(); 
  print OUT <<EOF;
</RedstarNPt>
EOF

  close OUT;

  return $file;  
}


#
#
#
# COLORVEC
#
#
#

sub print_colorvec_header
{
  &err_func(__func__,__LINE__,"expected 1 arg") unless $#_ == 1; 
  my ($self,$derivP) = @_; 

  print "starting : " . __func__ . "\n";

  my $t = $self->t_origin(); 
  my $n = $self->num_vecs(); 
  my $L = "false";
  my $S = "false";
  $L = "true" unless ! $self->convertUDtoL();
  $S = "true" unless ! $self->convertUDtoS();

  print OUT<<EOF;
<?xml version="1.0"?>
<ColorVecHadron>
  <Param>
    <version>4</version>
    <num_vecs>$n</num_vecs>
    <t_origin>$t</t_origin>
    <convertUDtoL>$L</convertUDtoL>
    <convertUDtoS>$S</convertUDtoS>
    <rcut>0.0</rcut>
    <use_derivP>$derivP</use_derivP>
EOF

  $self->print_flavor_to_mass_xml__generic();

  print OUT<<EOF;
  </Param>
EOF
  print "exiting : " . __func__ . "\n";
}

sub print_colorvec_db_files
{
  &err_func(__func__,__LINE__,"expected 7 args") unless $#_ == 7;
  my ($self,$node,$prop,$bar,$mes,$glue,$tet,$out) = @_;

  print OUT<<EOF;
  <DBFiles>
    <hadron_node_xmls>
EOF
  $self->print_elem_arr($node);
  print OUT<<EOF; 
    </hadron_node_xmls>
    <prop_dbs>
EOF
  $self->print_elem_arr($prop);
  print OUT<<EOF;
    </prop_dbs>
    <baryon_dbs>
EOF
  $self->print_elem_arr($bar); 
  print OUT<<EOF;
    </baryon_dbs>
    <meson_dbs>
EOF
  $self->print_elem_arr($mes);
  print OUT<<EOF;
    </meson_dbs>
    <glue_dbs>
EOF
  $self->print_elem_arr($glue);
  print OUT<<EOF;
    </glue_dbs>
    <tetra_dbs>
EOF
  $self->print_elem_arr($tet);
  print OUT<<EOF;
    </tetra_dbs>
    <output_db>$out</output_db>
  </DBFiles>
EOF

}

sub print_colorvec_trailer
{
  my $self = shift; 
  print OUT<<EOF;
</ColorVecHadron>
EOF
}

sub write_meson_hadron_node_ini_xml
{
  my $self = shift; 
  print __func__ . ": entering \n" unless ! $self->verbose();
  my @empty = (); 
  my $file = $self->scratch_seq_callback("meson_hadron_node"); 
  my $shadnode = $self->scratch_seq_callback("smeared_hadron_node_xml"); 
  my $shadop = $self->scratch_seq_callback("smeared_hadron_node_sdb");

  my @had_node_xmls = (); 
  push @had_node_xmls , $shadnode; 
  my $peram = $self->peram_dbs(); 
  my $mes = $self->meson_dbs(); 
  my $bar = $self->baryon_dbs(); 
  my $glue = $self->glue_dbs(); 
  my $tet = $self->tetraquark_dbs(); 
  open OUT , ">" , $file; 

  $self->print_colorvec_header("false"); 
  $self->print_colorvec_db_files(\@had_node_xmls,$peram,$bar,$mes,$glue,$tet,$shadop); 
  $self->print_colorvec_trailer();

  close OUT; 
  print __func__ . ": exiting \n" unless ! $self->verbose();
  return $file;
}

1;
