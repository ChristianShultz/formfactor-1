#!/usr/bin/perl

# make a params package since perl w/o using strict is scary 
package Params;
use File::Copy;
use File::Basename; 
use strict;


sub new
{
  my $invocant = shift;
  my $class = ref($invocant) || $invocant; 
  my $self =  {
    STEM => "szscl21_20_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per",
    BASEDIR => undef, 
    VERSION => 2,
    DIAGNOSTIC => 10,
    SEQNO => undef, 
    NUMVECS => 128,
    SNUMDILS => 128,
    TNUMDILS => 128, 
    NTCORR => 40,
    TSOURCE => [0],
    DT => [6],
    PROJOPXML => [],
    L_T =>128,
    L_S => 20,
    DECAYDIR => 3,
    CONVERTUDTOL => "true", 
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
    HADRONNODEXML => "hadron_node.xml",
    REDSTARXML=> "redstar.ini.xml",
    MESONHADRONNODE => "meson_hadron_node.ini.xml",
    HAROMHADRONNODE => "harom_hadron_node.ini.xml",
    GRAPHSDB => "graph.sdb",
    HADRONNODESDB => "hadron_node.sdb",
    OUTPUTSDB => "output.npt.sdb", 
    DISTSRCINI => "dist_src.ini.xml",
    PROPDISTINI => "prop_dist.ini.xml",
    PERAMDISTINI => "peram_dist.ini.xml",
    NPOINTLIST => "npt.list.xml",
    CONNQUARKLINEPOOL => [1,2,3,4,5,6,7,8,9,10],
    ANNIHQUARKLINEPOOL => [1,2,3,4,5,6,7,8],
    CLOVERCS => 1.589327,
    CLOVERCT => 0.902784,
    XI => 3.5,
    XI_0 => 4.3,
    NU => 1.265,
    RHO => 0.14,
    N_SMEAR => 2, 
    CFGTYPE => "SCIDAC",
    PROPDISTINVERTER => "BICGSTAB", 
    SKELETON=>undef,
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

sub version{
  my $self = shift; 
  if(@_) {$self->{VERSION} = shift}
  return $self->{VERSION};
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

sub num_space_dils{
  my $self = shift; 
  if(@_) {$self->{SNUMDILS} = shift}
  return $self->{SNUMDILS};
}

sub num_time_dils{
  my $self = shift; 
  if(@_) {$self->{TNUMDILS} = shift}
  return $self->{TNUMDILS};
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

sub hadron_node_xmls
{
  my $self = shift;
  if(@_) {$self->{HADRONNODEXML} = shift}
  return $self->{HADRONNODEXML};
}

sub redstar_ini {
  my $self = shift; 
  if(@_) {$self->{REDSTARXML} = shift}
  return $self->{REDSTARXML};
}

sub meson_hadron_node {
  my $self = shift; 
  if(@_) {$self->{MESONHADRONNODE} = shift}
  return $self->{MSEONHADRONNODE};
}

sub harom_hadron_node {
  my $self = shift; 
  if(@_) {$self->{HAROMHADRONNODE} = shift}
  return $self->{HAROMHADRONNODE};
}

sub graph_sdb {
  my $self = shift; 
  if(@_) {$self->{GRAPHSDB} = shift}
  return $self->{GRAPHSDB};
}

sub hadron_node_sdb {
  my $self = shift; 
  if(@_) {$self->{HADRONNODESDB} = shift}
  return $self->{HADRONNODESDB};
}

sub output_sdb {
  my $self = shift; 
  if(@_) {$self->{OUTPUTSDB} = shift}
  return $self->{OUTPUTSDB};
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

sub conn_quark_line_pool{
  my $self = shift; 
  if(@_) {$self->{CONNQUARKLINEPOOL} = shift}
  return $self->{CONNQUARKLINEPOOL};
}

sub annih_quark_line_pool{
  my $self = shift; 
  if(@_) {$self->{ANNIHQUARKLINEPOOL} = shift}
  return $self->{ANNIHQUARKLINEPOOL};
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

sub prop_dist_inverter{
  my $self = shift; 
  if(@_) {$self->{PROPDISTINVERTER} = shift}
  return $self->{PROPDISTINVERTER};
}

sub skeleton{
  my $self = shift; 
  if(@_) {$self->{SKELETON} = shift}
  return $self->{SKELETON};
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

sub prop_dbs{
  my $self = shift; 
  
  if ( $self->scratchy() )
  {
    
  if ( $self->verbose() ) 
  {
    my @t = @{$self->scratchy_props()};
    print "prop_dbs props @t \n";
  }

    return $self->scratchy_props(); 
  } 
  
  
  my @prop_dbs; 
  my $cache_dir = $self->cache_dir();
  my $stem = $self->stem();
  my $seqno = $self->seqno();
  my $t_ref = $self->t_sources(); 
  my @flav = ("light");
  foreach my $tt (@$t_ref)
  {
    foreach my $flavor (@flav)
    {
      push @prop_dbs , "$cache_dir/${stem}/prop_db/${stem}.prop.n128.${flavor}.t0_${tt}.sdb${seqno}";
    }
  }

# don't need diag guys since no bubbles

  my @scratch_dbs = ();

  if($self->skeleton())
  {
    @scratch_dbs = @{ $self->gen_prop_peram_files("all") };
  }
  else
  {
    foreach my $db (@prop_dbs)
    {
      push @scratch_dbs , $self->copy_file_to_scratch($db);
    }
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
  
  if ( $self->scratchy() ) 
  {
    return $self->scratchy_perams(); 
  }
  else
  {
    print " I DONT KNOW WHAT TO DO \n ";
    die ( " \nBOOM \n\n" );
  }
  my @stub = (); 
  return \@stub;
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

  my $lhpc_dir = $self->lhpc_dir(); 
  my $lhpc_stem = "$lhpc_dir/${stem}/meson_db_128/${stem}.meson.colorvec.sdb";
  push @mdbs , $lhpc_stem . $self->seqno(); 

  my @scratch_dbs = ();

  if($self->skeleton())
  {
    push @scratch_dbs , "some_database" ; 
  }
  else
  {
    foreach my $db (@mdbs)
    {
      push @scratch_dbs , $self->copy_file_to_scratch($db); 
    }
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

#
#
#  call back functions for filenames that know about skeleton-ness
#
#


sub sequenceify 
{
  my $self = shift; 
  my $name = shift; 

  if ($self->skeleton())
  {
    $name = $name . "1";
  }
  else
  {
    $name = $name . $self->seqno(); 
  }

  return $name; 
}

sub get_graph_db_string
{
  my $self = shift; 
  my $x = $self->graph_sdb(); 
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_hadron_node_xml_string
{
  my $self = shift; 
  my $x = $self->hadron_node_xmls(); 
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x);
}

sub get_redstar_xml_ini_string
{
  my $self = shift;
  my $x = $self->redstar_ini();
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_meson_hadron_node_xml_ini_string
{
  my $self = shift; 
  my $x = $self->meson_hadron_node();
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_harom_hadron_node_xml_ini_string
{
  my $self = shift;
  my $x = $self->harom_hadron_node(); 
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_hadron_node_output_db_string
{
  my $self = shift; 
  my $x = $self->hadron_node_sdb();
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_output_db_string
{
  my $self = shift; 
  my $x = $self->output_sdb();
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x);  
}

sub get_dist_src_ini_string
{
  my $self = shift; 
  my $x = $self->dist_src_ini(); 
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_prop_dist_ini_string
{
  my $self = shift; 
  my $x = $self->prop_dist_ini(); 
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub get_peram_dist_ini_string
{
  my $self = shift; 
  my $x = $self->peram_dist_ini();
  $x = $self->scratchify_file_name($x);  
  return $self->sequenceify($x); 
}

sub proj_ops_xml
{
  my $self = shift; 
  my $ref = $self->{PROJOPXML};

  if($self->sequenceify())
  {
    my @foo = (); 
    return \@foo; 
  }

  return $ref; 
}

#
#
#  PRINITERS
#
#









#
#
#
#  COLORVEC STUFF
#
#
#

## DEPRECIATED -- not in use for unsmeared gen props

sub print_meson_hadron_node_Param
{
  my $self = shift; 
  my $t_origin = $self->t_origin();

  my $um = $self->mass_string($self->u_mass()); 
  my $sm = $self->mass_string($self->s_mass()); 

  print OUT <<EOF;
  <Param>
    <version>$self->{VERSION}</version>
    <num_vecs>$self->{NUMVECS}</num_vecs>
    <Nt_corr>$self->{NTCORR}</Nt_corr>
    <t_origin>$t_origin</t_origin>

    <FlavorToMass>
      <elem>
        <flavor>l</flavor>
        <mass>$um</mass>
      </elem>
      <elem>
        <flavor>s</flavor>
        <mass>$sm</mass>
      </elem>
    </FlavorToMass>
  </Param>
EOF
}

sub print_meson_hadron_node_DBFiles
{
  my $self = shift; 
  print OUT <<EOF;
  <DBFiles>
    <hadron_node_xmls>
      <elem>/scratch/$self->{HADRONNODEXML}</elem>
    </hadron_node_xmls>
    <prop_dbs>
EOF
  my $ref = $self->prop_dbs(); 
  my @prop_dbs = @$ref;
  foreach (@prop_dbs)
  {
    print OUT <<EOF;
      <elem>$_</elem>
EOF
  }

  print OUT <<EOF;
    </prop_dbs>

    <baryon_dbs>
EOF
  my $ref = $self->baryon_dbs();
  my @dbs = @$ref;
  foreach (@dbs)
  {
    print OUT <<EOF;
      <elem>$_</elem>
EOF
  }

  print OUT <<EOF;
    </baryon_dbs>

    <meson_dbs>
EOF
  my $ref = $self->meson_dbs();
  my @dbs = @$ref;
  foreach (@dbs)
  {
    print OUT <<EOF;
      <elem>$_</elem>
EOF
  }

  print OUT <<EOF;
    </meson_dbs>

    <glue_dbs>
EOF
  my $ref = $self->glue_dbs();
  my @dbs = @$ref;
  foreach (@dbs)
  {
    print OUT <<EOF;
      <elem>$_</elem>
EOF
  }

  my $out_db = $self->hadron_node_output_db(); 
  print OUT <<EOF;
    </glue_dbs>

    <output_db>$out_db</output_db>
  </DBFiles>
EOF
}

## DONT USE ME !!!! 

sub write_meson_hadron_node_xml
{
  my $self = shift; 
  my $file = $self->get_meson_hadron_node_xml_ini_string();


  if(-e $file)
  {
    unlink($file) || die($_); 
  }

  open OUT , ">" , $file; 

  print OUT <<EOF;
<?xml version="1.0"?>
<ColorVecHadron>
EOF
  $self->print_meson_hadron_node_Param(); 
  $self->print_meson_hadron_node_DBFiles();  
  print OUT <<EOF;
</ColorVecHadron>
EOF
  close OUT; 

  return $file; 
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
    <version>$self->{VERSION}</version>
    <diagnostic_level>$self->{DIAGNOSTIC}</diagnostic_level>
    <Nt_corr>$self->{NTCORR}</Nt_corr>
    <convertUDtoL>$self->{CONVERTUDTOL}</convertUDtoL>
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
  my $proj_ops = $self->proj_ops_xml(); 
  my @projops = @$proj_ops; 
  my $had_graph = $self->get_graph_db_string(); 
  my $had_node_xml = $self->get_hadron_node_xml_string(); 
  my $had_node_out_db = $self->get_hadron_node_output_db_string();
  my $out_db = $self->get_output_db_string();  
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
    <hadron_node_xml>$had_node_xml</hadron_node_xml>
    <hadron_npt_graph_db>$had_graph</hadron_npt_graph_db>
    <hadron_node_dbs><elem>$had_node_out_db</elem></hadron_node_dbs>
    <output_db>$out_db</output_db>
  </DBFiles>
EOF
}

sub write_redstar_xml
{
  my $self = shift;

  my $file = $self->get_redstar_xml_ini_string(); 

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

sub colorvec_file
{
  my $self = shift; 

  my $c = ""; 

  if($self->skeleton())
  {
    $c = "./colorvec.timeslice.mod"; 
  }
  else
  {
    my $cf = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->stem() . "/eigs_mod/" . $self->stem() . ".3d.eigs.mod" . $self->seqno(); 
    $c = $self->copy_file_to_scratch($cf); 
  }

  return $c; 
}


sub gauge_file_lime
{
  my $self = shift; 
  my $g = ""; 
  if( $self->skeleton() )
  {
    return $self->gauge_file_3d(); 
  }
  else
  { 
    my $path = "/lustre/cache/LHPC/Spectrum/Clover/NF2+1/" . $self->stem() . "/cfgs/";
    my $gf = $path . $self->stem() . "_cfg_" . $self->seqno() . ".lime"; 
    $g = $self->copy_file_to_scratch($gf); 
  }

  return $g;
}

sub gauge_file_3d
{
  my $self = shift; 
  my $g = "";

  if($self->skeleton())
  {
    $g = "./gauge.timeslice.mod"; 
  }  
  else
  {
    my $gfpath = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->stem() . "/cfgs_mod/";
    my $gf = $gfpath . $self->stem() . ".3d.gauge.mod" . $self->seqno(); 
    $g = $self->copy_file_to_scratch($gf); 
  }
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

sub quark_line_string
{
  my ($self, $delim , $line_type) = @_;   

  if($line_type eq "CONN")
  {
    return $self->arr_string($delim,$self->conn_quark_line_pool());
  }
  elsif($line_type eq "ANNIH")
  {
    return $self->arr_string($delim,$self->annih_quark_line_pool()); 
  }
  else
  {
    print "ERROR: unknown quark line type, $line_type , exiting.. \n";
    die("unholy death");  
  }

  die($_); 
}

sub time_source_string
{
  my $self = shift;
  my $delim = shift; 
  my $line_type = shift;  


  print " time_source_string delim = 1 $delim 1 , line = $line_type \n" unless ! $self->verbose();

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


sub gen_prop_src_files
{
  my $self = shift;
  my $opt = shift;
  my $line_type = shift; 


  print "gen_prop_src_files opt = $opt , line = $line_type \n";

  my @tmp = (); 
  my @ret = (); 

  my  $quark_line_string = ".q" . $self->quark_line_string("_", $line_type); 

  if($self->skeleton())
  {
    if ($opt eq "f")
    {
      push @tmp , "srcf.mod"; 
    }
    elsif ( $opt eq "b" )
    {
      push @tmp , "srcb.mod"; 
    }
    elsif ( $opt eq "a" )
    {
      push @tmp , "srca.mod"; 
    }
    elsif ($opt eq "all")
    {
      push @tmp , "srcf.mod"; 
      push @tmp , "srcb.mod"; 
      push @tmp , "srca.mod"; 
    }
    else
    {
      print "ERROR: unrecognized option $opt, exiting.. \n";
      die("unholy death");
    }
  }
  else
  {
    my $stemy = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->ensemble() . "/gen_props/soln_vectors";
    my $tail = "." . $line_type . ".src.n" . $self->num_vecs() . ".d" . $self->num_space_dils() ;

    if ($opt eq "f")
    {    
      my $f = $stemy . "/forward/";
      $tail .= $quark_line_string . ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $f . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "b")
    {
      my $b = $stemy . "/backward/";
      $tail .= $quark_line_string . ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $b . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "a")
    {
      my $b = $stemy . "/annih/";
      $tail .= $quark_line_string . ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $b . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "all")
    {
      if ($line_type != "CONN") 
      {
        die ("NEED TO USE CONN HERE \n"); 
      }
      my  $quark_line_string_a = ".q" . $self->quark_line_string("_", "ANNIH"); 
      my $b = $stemy . "/backward/";
      push @tmp , $b . $self->ensemble() . $tail . $quark_line_string . ".t0" . $self->time_source_string("_","b") . ".mod"; 
      my $f = $stemy . "/forward/";
      push @tmp , $f . $self->ensemble() . $tail . $quark_line_string . ".t0" . $self->time_source_string("_","f") . ".mod";
      my $a = $stemy . "/annih/";
      push @tmp , $a . $self->ensemble() . $tail . $quark_line_string_a . ".t0" . $self->time_source_string("_","a") . ".mod";
    }
    else
    {
      print "ERROR: unrecognized option $opt, exiting.. \n";
      die("unholy death");
    }

  }

  foreach my $elem (@tmp)
  {
    push @ret, $self->sequenceify($elem)
  }

  return \@ret; 
}



sub gen_prop_soln_vector_files
{
  my $self = shift;
  my $opt = shift;
  my $line_type = shift; 

  print "gen_prop_soln_vector_files opt = $opt line_type = $line_type \n";

  my @tmp = (); 
  my @ret = (); 

  my  $quark_line_string = ".q" . $self->quark_line_string("_", $line_type); 


  if($self->skeleton())
  {
    if ($opt eq "f")
    {
      push @tmp , "solnf.mod"; 
    }
    elsif ( $opt eq "b" )
    {
      push @tmp , "solnb.mod"; 
    }
    elsif ($opt eq "a" )
    {
      # NOTE HACK
      push @tmp , "solnf.mod";
    }
    elsif ($opt eq "all")
    {
      push @tmp , "solnf.mod"; 
      push @tmp , "solnb.mod"; 
    }
    else
    {
      print "ERROR: unrecognized option $opt, exiting.. \n";
      die("unholy death");
    }
  }
  else
  {
    my $stemy = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->ensemble() . "/gen_props/soln_vectors";
    my $tail = "." . $line_type . ".prop.n" . $self->num_vecs() . ".d" . $self->num_space_dils() ;

    if ($opt eq "f")
    {    
      my $f = $stemy . "/forward/";
      $tail .= $quark_line_string . ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $f . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "b")
    {
      my $b = $stemy . "/backward/";
      $tail .= $quark_line_string . ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $b . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "a")
    {
      # NOTE HACK
      my $f = $stemy . "/forward/";
      $tail .= $quark_line_string . ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $f . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "all")
    {
      if ($line_type != "CONN") 
      {
        die ("NEED TO USE CONN HERE \n"); 
      }
      my $b = $stemy . "/backward/";
      push @tmp , $b . $self->ensemble() . $tail . $quark_line_string . ".t0" . $self->time_source_string("_","b") . ".mod"; 
      my $f = $stemy . "/forward/";
      push @tmp , $f . $self->ensemble() . $tail . $quark_line_string . ".t0" . $self->time_source_string("_","f") . ".mod";
    }
    else
    {
      print "ERROR: unrecognized option $opt, exiting.. \n";
      die("unholy death");
    }

  }

  foreach my $elem (@tmp)
  {
    push @ret, $self->sequenceify($elem)
  }

  return \@ret; 
}


sub gen_prop_peram_files
{
  my $self = shift;
  my $opt = shift;
  my $line_type = shift; 

  print "gen_prop_peram_files opt = $opt , line = $line_type " unless ! $self->verbose();

  my @tmp = (); 
  my @ret = (); 


  if($self->skeleton())
  {
    if ($opt eq "f")
    {
      push @tmp , "peramf.mod"; 
    }
    elsif ( $opt eq "b" )
    {
      push @tmp , "peramb.mod"; 
    }
    elsif ( $opt eq "a" )
    {
      push @tmp , "perama.mod"; 
    }
    elsif ($opt eq "all")
    {
      push @tmp , "peramf.mod"; 
      push @tmp , "peramb.mod";
      push @tmp , "perama.mod"; 
    }
    else
    {
      print "ERROR: unrecognized option $opt, exiting.. \n";
      die("unholy death");
    }
  }
  else
  {
    my $stemy = "/lustre/cache/Spectrum/Clover/NF2+1/" . $self->ensemble() . "/gen_props/soln_vectors";
    my $tail = "." . $line_type . ".peram.n" . $self->num_vecs() . ".d" . $self->num_space_dils() ;

    if ($opt eq "f")
    {    
      my $f = $stemy . "/forward/";
      $tail .= ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $f . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "b")
    {
      my $b = $stemy . "/backward/";
      $tail .=  ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $b . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "a")
    {
      my $b = $stemy . "/annih/";
      $tail .=  ".t0" . $self->time_source_string("_",$opt) . ".mod";   
      push @tmp , $b . $self->ensemble() . $tail; 
    }
    elsif ($opt eq "all")
    {
      if ($line_type != "CONN") 
      {
        die ("NEED TO USE CONN HERE \n"); 
      }
      my $b = $stemy . "/backward/";
      push @tmp , $b . $self->ensemble() . $tail . ".t0" . $self->time_source_string("_","b") . ".mod"; 
      my $f = $stemy . "/forward/";
      push @tmp , $f . $self->ensemble() . $tail . ".t0" . $self->time_source_string("_","f") . ".mod";
      my $a = $stemy . "/annih/";
      my $tail_a = "." . "ANNIH" . ".peram.n" . $self->num_vecs() . ".d" . $self->num_space_dils() ;
      push @tmp , $a . $self->ensemble() . $tail_a . ".t0" . $self->time_source_string("_","a") . ".mod";
    }
    else
    {
      print "ERROR: unrecognized option $opt, exiting.. \n";
      die("unholy death");
    }

  }

  foreach my $elem (@tmp)
  {
    push @ret, $self->sequenceify($elem)
  }

  return \@ret; 
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


sub print_distillution_noise
{
  my $self = shift; 

  print OUT <<EOF; 

      <elem>
        <annotation>Distillution noise</annotation>
        <Name>DISTILLUTION_NOISE</Name>
        <Frequency>1</Frequency>
        <Param>
          <ensemble>$self->{STEM}</ensemble>
          <sequence>$self->{SEQNO}</sequence>
          <decay_dir>$self->{DECAYDIR}</decay_dir>
        </Param>
        <NamedObject>
          <distillution_id>dist_obj</distillution_id>
        </NamedObject>
      </elem>

EOF
}

sub print_quark_line_xml
{
  my ($self,$line_type,$tsources,$ntf,$ntb) = @_; 

  if ( $line_type eq "CONN" )
  {
    print OUT <<EOF;
          <QuarkLine>
           <QuarkLineType>CONN</QuarkLineType>
            <num_vecs>$self->{NUMVECS}</num_vecs>
            <num_space_dils>$self->{SNUMDILS}</num_space_dils>
            <t_sources>$tsources</t_sources>
            <Nt_forward>$ntf</Nt_forward>
            <Nt_backward>$ntb</Nt_backward>
          </QuarkLine>
EOF
  }
  elsif ($line_type eq "ANNIH" )
  {
    print OUT <<EOF;
          <QuarkLine>
            <QuarkLineType>ANNIH</QuarkLineType>
            <num_vecs>$self->{NUMVECS}</num_vecs>
            <num_space_dils>$self->{SNUMDILS}</num_space_dils>
            <num_time_dils>$self->{TNUMDILS}</num_time_dils>
          </QuarkLine>
EOF
  }
  else
  {
    print "ERROR: unknown line type $line_type, exiting.. ";
    die("unholy_death");
  }


}


sub print_make_source_dist
{
  my ($self , $mass, $line_type , $f_or_b ) = @_; 

  print " print_make_source_dist mass = $mass , line = $line_type , forb = $f_or_b \n" unless ! $self->verbose(); 

  $mass = $self->mass_string($mass); 

  my $lines = substr( $self->quark_line_string(" ",$line_type) , 1);   #beautification
  my $tsources = substr( $self->time_source_string(" ", $f_or_b ) , 1);  
  my $ntf = $self->nt_corr(); 
  my $ntb = $self->nt_corr(); 

  if($f_or_b eq "f")
  {
    $ntb = 0; 
  }
  elsif ($f_or_b eq "b") 
  {
    $ntf = 0; 
  }
  elsif ($f_or_b eq "a")
  {
    $ntb = 0; 
  }
  else
  {
    die ("unrecognized opt , $f_or_b \n"); 
  }

  my $colorvec = $self->colorvec_file(); 
  my $src_ref = $self->gen_prop_src_files($f_or_b, $line_type);

  my $srcfull = ${$src_ref}[0]; 
  my @src_files = @{$self->scratchify_file_name_arr($src_ref) };
  die("too many source files ( $#src_files ) " ) unless ($#src_files == 0); 
  my $src = $src_files[0]; 


  print " tsources = $tsources , lines = $lines , src = $src " unless ! $self->verbose(); 

  print OUT <<EOF; 

    <elem>
      <Name>MAKE_SOURCE_DIST</Name>
      <Frequency>1</Frequency>
      <Param>
      <Contractions>
        <quark_lines>$lines</quark_lines>
        <mass>$mass</mass>
EOF

  $self->print_quark_line_xml($line_type,$tsources,$ntf,$ntb); 

  print OUT <<EOF; 
      </Contractions>
    </Param>
    <NamedObject>
      <distillution_id>dist_obj</distillution_id>
      <colorvec_file>$colorvec</colorvec_file>
      <src_file>$src</src_file>
    </NamedObject>
  </elem>

EOF

  return $srcfull; 
}

#
sub print_chroma_prop_distillution_BICGSTAB
{
  my ($self , $mass , $line_type , $f_or_b) = @_; 

  my $m_q = $mass; 

  $mass = $self->mass_string($mass);

  my $lines = substr( $self->quark_line_string(" ",$line_type) , 1);   #beautification
  my $tsources = substr( $self->time_source_string(" ", $f_or_b ) , 1);  
  my $ntf = $self->nt_corr(); 
  my $ntb = $self->nt_corr(); 

  if($f_or_b eq "f")
  {
    $ntb = 0; 
  }
  elsif ($f_or_b eq "b") 
  {
    $ntf = 0; 
  }
  elsif ($f_or_b eq "a")
  {
    $ntb = 0; 
  }
  else
  {
    die ("unrecognized opt , $f_or_b \n"); 
  }

  my $colorvec = $self->colorvec_file(); 
  my $src_ref = $self->gen_prop_src_files($f_or_b, $line_type);

  my @src_files = @{$self->scratchify_file_name_arr($src_ref) };
  die("too many source files ( $#src_files ) " ) unless ($#src_files == 0); 
  my $src = $src_files[0]; 

  my $prop_ref = $self->gen_prop_soln_vector_files($f_or_b , $line_type); 
  my $propfull = ${$prop_ref}[0]; 
  my @prop_files = @{$self->scratchify_file_name_arr($prop_ref) };
  die("too many prop files ( $#prop_files ) " ) unless ($#prop_files == 0); 
  my $prop = $prop_files[0];


  print OUT <<EOF;
  <elem>
    <Name>PROP_DISTILLUTION</Name>
    <Frequency>1</Frequency>
    <Param>
      <Contractions>
        <quark_lines>$lines</quark_lines>
        <mass>$mass</mass>
EOF

  $self->print_quark_line_xml($line_type,$tsources,$ntf,$ntb); 


  print OUT <<EOF;
        </Contractions>
        <Propagator>
          <version>10</version>
          <quarkSpinType>FULL</quarkSpinType>
          <obsvP>false</obsvP>
          <numRetries>1</numRetries>
          <FermionAction>
            <FermAct>CLOVER</FermAct>
            <Mass>$m_q</Mass>
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
          <InvertParam>
            <invType>BICGSTAB_INVERTER</invType>
            <RsdBiCGStab>1.e-7</RsdBiCGStab>
            <MaxBiCGStab>1000</MaxBiCGStab>
          </InvertParam>
         </Propagator>
      </Param>
      <NamedObject>
        <save_solnP>true</save_solnP>
        <gauge_id>default_gauge_field</gauge_id>
        <distillution_id>dist_obj</distillution_id>
        <src_file>$src</src_file>
        <soln_file>$prop</soln_file>
      </NamedObject>
    </elem>

EOF

  return $propfull; 
}


sub print_chroma_prop_distillution_CG_TEST
{
  my ($self , $mass , $line_type , $f_or_b) = @_; 

  my $m_q = $mass;

  $mass = $self->mass_string($mass); 

  print " print_chroma_prop_distillution mq = $mass , line = $line_type , opt = $f_or_b \n";

  my $lines = substr( $self->quark_line_string(" ",$line_type) , 1);   #beautification
  my $tsources = substr( $self->time_source_string(" ", $f_or_b ) , 1);  
  my $ntf = $self->nt_corr(); 
  my $ntb = $self->nt_corr(); 


  if($f_or_b eq "f")
  {
    $ntb = 0; 
  }
  elsif ($f_or_b eq "b") 
  {
    $ntf = 0; 
  }
  elsif ($f_or_b eq "a")
  {
    $ntb = 0; 
  }
  else
  {
    die ("unrecognized opt , $f_or_b \n"); 
  }

  my $colorvec = $self->colorvec_file(); 
  my $src_ref = $self->gen_prop_src_files($f_or_b, $line_type);

  my @src_files = @{$self->scratchify_file_name_arr($src_ref) };
  die("too many source files ( $#src_files ) " ) unless ($#src_files == 0); 
  my $src = $src_files[0]; 

  my $prop_ref = $self->gen_prop_soln_vector_files($f_or_b , $line_type); 
  my $propfull = ${$prop_ref}[0];
  my @prop_files = @{$self->scratchify_file_name_arr($prop_ref) };
  die("too many prop files ( $#prop_files ) " ) unless ($#prop_files == 0); 
  my $prop = $prop_files[0];


  print OUT <<EOF;
  <elem>
    <Name>PROP_DISTILLUTION</Name>
    <Frequency>1</Frequency>
    <Param>
      <Contractions>
        <quark_lines>$lines</quark_lines>
        <mass>$mass</mass>
EOF

  $self->print_quark_line_xml($line_type,$tsources,$ntf,$ntb); 


  print OUT <<EOF;
        </Contractions>
        <Propagator>
          <version>10</version>
          <quarkSpinType>FULL</quarkSpinType>
          <obsvP>false</obsvP>
          <numRetries>1</numRetries>
          <FermionAction>
            <FermAct>CLOVER</FermAct>
            <Mass>$m_q</Mass>
            <clovCoeff>1.0</clovCoeff>
            <AnisoParam>
              <anisoP>false</anisoP>
            </AnisoParam>
            <FermionBC>
              <FermBC>SIMPLE_FERMBC</FermBC>
              <boundary>1 1 1 -1</boundary>
            </FermionBC>
          </FermionAction>
          <InvertParam>
            <invType>CG_INVERTER</invType>
            <RsdCG>1.e-8</RsdCG>
            <MaxCG>1000</MaxCG>
          </InvertParam>
         </Propagator>
      </Param>
      <NamedObject>
        <save_solnP>true</save_solnP>
        <gauge_id>default_gauge_field</gauge_id>
        <distillution_id>dist_obj</distillution_id>
        <src_file>$src</src_file>
        <soln_file>$prop</soln_file>
      </NamedObject>
    </elem>

EOF

  return $propfull; 
}



sub print_chroma_prop_distillution_QUDA
{
  my ($self , $mass , $line_type , $f_or_b) = @_; 

  my $m_q = $mass; 

  $mass = $self->mass_string($mass);

  print " print_chroma_prop_distillution mq = $mass , line = $line_type , opt = $f_or_b \n";

  my $lines = substr( $self->quark_line_string(" ",$line_type) , 1);   #beautification
  my $tsources = substr( $self->time_source_string(" ", $f_or_b ) , 1);  
  my $ntf = $self->nt_corr(); 
  my $ntb = $self->nt_corr(); 

  if($f_or_b eq "f")
  {
    $ntb = 0; 
  }
  elsif ($f_or_b eq "b") 
  {
    $ntf = 0; 
  }
  elsif ($f_or_b eq "a")
  {
    $ntb = 0; 
  }
  else
  {
    die ("unrecognized opt , $f_or_b \n"); 
  }

  my $colorvec = $self->colorvec_file(); 
  my $src_ref = $self->gen_prop_src_files($f_or_b, $line_type);

  my @src_files = @{$self->scratchify_file_name_arr($src_ref) };
  die("too many source files ( $#src_files ) " ) unless ($#src_files == 0); 
  my $src = $src_files[0]; 

  my $prop_ref = $self->gen_prop_soln_vector_files($f_or_b , $line_type); 
  my $propfull = ${$prop_ref}[0];
  my @prop_files = @{$self->scratchify_file_name_arr($prop_ref) };
  die("too many prop files ( $#prop_files ) " ) unless ($#prop_files == 0); 
  my $prop = $prop_files[0];

  print OUT <<EOF;
  <elem>
    <Name>PROP_DISTILLUTION</Name>
    <Frequency>1</Frequency>
    <Param>
      <Contractions>
        <quark_lines>$lines</quark_lines>
        <mass>$mass</mass>
EOF

  $self->print_quark_line_xml($line_type,$tsources,$ntf,$ntb); 

  print OUT <<EOF;
      </Contractions>
      <Propagator>
        <version>10</version>
        <quarkSpinType>FULL</quarkSpinType>
        <obsvP>false</obsvP>
        <numRetries>1</numRetries>
        <FermionAction>
          <FermAct>CLOVER</FermAct>
          <Mass>$m_q</Mass>
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
        <InvertParam>
          <invType>QUDA_CLOVER_INVERTER</invType>
          <CloverParams>
            <Mass>$m_q</Mass>
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
      </Propagator>
    </Param>
    <NamedObject>
      <save_solnP>true</save_solnP>
      <gauge_id>default_gauge_field</gauge_id>
      <distillution_id>dist_obj</distillution_id>
      <src_file>$src</src_file>
      <soln_file>$prop</soln_file>
    </NamedObject>
  </elem>

EOF

  return $propfull; 
}


sub print_chroma_prop_distillution
{
  my ($self,$mass,$line_type,$f_or_b) = @_; 

  my $jango = $self->prop_dist_inverter(); 

  print " print_chroma_prop_distillution mass = $mass , line = $line_type ,";
  print " forb = $f_or_b , inv = $jango \n";

  my %hash = (); 
  $hash{"BICG_STAB"} = \&Params::print_chroma_prop_distillution_BICGSTAB; 
  $hash{"CG_TEST"} = \&Params::print_chroma_prop_distillution_CG_TEST; 
  $hash{"QUDA"} = \&Params::print_chroma_prop_distillution_QUDA;

  if ( exists $hash{$self->prop_dist_inverter()} )
  {
    return  &{ $hash{$self->prop_dist_inverter()} }($self,$mass,$line_type,$f_or_b);
  }
  else
  {
    my @k = keys %hash; 
    print "ERROR: could not find " . $self->prop_dist_inverter() . " in hash, choices are: \n";

    foreach my $kk (@k) 
    {
      print $kk . "\n"; 
    }
    die("unholy death"); 
  }
}



sub print_harom_peram_distillution
{

  my ($self , $mass , $line_type , $f_or_b) = @_; 

  $mass = $self->mass_string($mass);

  print " print_haroma_peram_distillution mq = $mass , line = $line_type , opt = $f_or_b \n";

  my $lines = substr( $self->quark_line_string(" ",$line_type) , 1);   #beautification
  my $tsources = substr( $self->time_source_string(" ", $f_or_b ) , 1);  
  my $ntf = $self->nt_corr(); 
  my $ntb = $self->nt_corr(); 

  if($f_or_b eq "f")
  {
    $ntb = 0; 
  }
  elsif ($f_or_b eq "b")
  {
    $ntf = 0; 
  }
  elsif ($f_or_b eq "a")
  {
    $ntb = 0; 
  }
  else
  {
    die("unrecognized opt : $f_or_b \n");
  }

  my $colorvec = $self->colorvec_file(); 
  my $src_ref = $self->gen_prop_src_files($f_or_b, $line_type);

  my @src_files = @{$self->scratchify_file_name_arr($src_ref) };
  die("too many source files ( $#src_files ) " ) unless ($#src_files == 0); 
  my $src = $src_files[0]; 

  my $prop_ref = $self->gen_prop_soln_vector_files($f_or_b , $line_type); 
  my @prop_files = @{$self->scratchify_file_name_arr($prop_ref) };
  die("too many prop files ( $#prop_files ) " ) unless ($#prop_files == 0); 
  my $prop = $prop_files[0];

  my $peram_ref = $self->gen_prop_peram_files($f_or_b , $line_type); 
  my $peramfull = ${$peram_ref}[0];
  my @peram_files = @{$self->scratchify_file_name_arr($peram_ref)};
  die("too many prop files ( $#peram_files ) " ) unless ($#peram_files == 0); 
  my $peram = $peram_files[0];

  print OUT <<EOF;
  <elem>
    <Name>PERAM_DIST</Name>
    <Frequency>1</Frequency>
    <Param>
      <Contractions>
        <quark_lines>$lines</quark_lines>
        <mass>$mass</mass>
EOF

  $self->print_quark_line_xml($line_type,$tsources,$ntf,$ntb); 

  print OUT <<EOF;
        </Contractions>
      </Param>
      <NamedObject>
        <distillution_id>dist_obj</distillution_id>
        <colorvec_file>$colorvec</colorvec_file>
        <soln_file>$prop</soln_file>
        <peram_file>$peram</peram_file>
      </NamedObject>
    </elem>

EOF

 return $peramfull; 
}



#
#
#  MAKE A DISTILLATION_SOURCE
#
#

sub write_dist_src_ini_xml
{
  my $self = shift; 
  my $file = $self->get_dist_src_ini_string(); 


  if (-e $file)
  { 
    unlink ($file) || die($_); 
  }

  open OUT , ">" , $file;

  $self->print_harom_header();

  $self->print_distillution_noise();

  my @ffs = (); 

  push @ffs, $self->print_make_source_dist($self->u_mass(),"CONN","f"); 
  push @ffs, $self->print_make_source_dist($self->u_mass(),"CONN","b");
  push @ffs, $self->print_make_source_dist($self->u_mass(),"ANNIH","a");

  $self->print_harom_trailer(); 

  close OUT;

  my @r ;
  push @r , $file , \@ffs; 

  
  if ( $self->verbose() ) 
  {
    print " dist_src.ini.xml = $file \n  all the outputs \n";
    foreach my $d (@ffs)
    {
      print $d . "\n";
    }
  }


  return \@r ;  
}

#
#
#
#
#  COMPUTE SOLUTION VECTORS WITH CHROMA 
#
#
#
#
#
#

sub write_prop_dist_chroma_ini_xml
{
  my $self = shift; 
  my $file = $self->get_prop_dist_ini_string(); 

  if (-e $file)
  {
    unlink ($file) || die($_); 
  }

  open OUT , ">" , $file ; 

  $self->print_chroma_header();

  $self->print_distillution_noise(); 

  my @ffs = (); 

  push @ffs, $self->print_chroma_prop_distillution($self->u_mass(),"CONN","f");
  push @ffs, $self->print_chroma_prop_distillution($self->u_mass(),"CONN","b"); 

  $self->print_chroma_trailer(); 

  close OUT ;

  my @r = (); 
  push @r, $file , \@ffs; 

  return \@r ; 
}

#
#
#
#
#  COMPUTE PERAMBULATORS WITH HAROM 
#
#
#


sub write_peram_dist_harom_ini_xml
{
  my $self = shift ; 
  my $file = $self->get_peram_dist_ini_string(); 

  if ( -e $file ) 
  {
    unlink $file || die ($_); 
  }

  open OUT , ">" , $file; 

  $self->print_harom_header(); 
  $self->print_distillution_noise(); 

  my @ffs = (); 
  
  push @ffs, $self->print_harom_peram_distillution($self->u_mass(),"CONN","f"); 
  push @ffs, $self->print_harom_peram_distillution($self->u_mass(),"CONN","b");
  push @ffs, $self->print_harom_peram_distillution($self->u_mass(),"ANNIH","a");

  $self->print_harom_trailer(); 

  close OUT;

  my @r = (); 
  push @r , $file , \@ffs; 

  return \@r;  
}


#
#
#
#
#  HAROM HADRON NODE 
#
#
#
#
sub print_harom_hadron_node_dist
{
  my $self = shift; 

  my $gauge = $self->gauge_file_3d(); 
  my $colorvec = $self->colorvec_file(); 
  my $prop_files = $self->prop_dbs();
  my $peram_files =  $self->peram_dbs() ;
  my $hadron_node_xml = $self->get_hadron_node_xml_string(); 
  my $hadron_node_out_db = $self->get_hadron_node_output_db_string();

  my $um = $self->mass_string($self->u_mass()); 
  my $sm = $self->mass_string($self->s_mass()); 
  my $conn_lines = $self->quark_line_string(" ","CONN"); 
  my $annih_lines = $self->quark_line_string(" ","ANNIH"); 


  if ( $self->verbose() ) 
  {
    my @t = @{$prop_files};
    print " props @t \n";
    @t = @{$peram_files}; 
    print " perams @t \n"; 
  }


  print OUT <<EOF; 

      <elem>
        <annotation>propagator soln vectors</annotation>
        <Name>HADRON_NODE</Name>
        <Frequency>1</Frequency>
        <Param>
          <version>1</version>
          <num_vecs>$self->{NUMVECS}</num_vecs>
          <conn_num_dils>$self->{SNUMDILS}</conn_num_dils>
          <annih_num_dils>$self->{SNUMDILS}</annih_num_dils>
          <displacement_length>1</displacement_length>

          <FlavorToMass>
            <elem>
              <flavor>l</flavor>
              <mass>$um</mass>
            </elem>
            <elem>
              <flavor>s</flavor>
              <mass>$sm</mass>
            </elem>
          </FlavorToMass>

          <ConnQuarkLinePool>
            <elem>
              <flavor>l</flavor>
              <quark_lines>$conn_lines</quark_lines>
            </elem>
            <elem>
              <flavor>s</flavor>
              <quark_lines>$conn_lines</quark_lines>
            </elem>
          </ConnQuarkLinePool>

          <AnnihQuarkLinePool>
            <elem>
              <flavor>l</flavor>
              <quark_lines>$annih_lines</quark_lines>
            </elem>
            <elem>
              <flavor>s</flavor>
              <quark_lines>$annih_lines</quark_lines>
            </elem>
          </AnnihQuarkLinePool>

          <LinkSmearing>
            <LinkSmearingType>$self->{LINKSMEARTYPE}</LinkSmearingType>
            <link_smear_fact>$self->{LINKSMEARFAC}</link_smear_fact>
            <link_smear_num>$self->{LINKSMEARNUM}</link_smear_num>
            <no_smear_dir>$self->{DECAYDIR}</no_smear_dir>
          </LinkSmearing>
        </Param>

        <NamedObject>
          <distillution_id>dist_obj</distillution_id>
          <gauge_file>$gauge</gauge_file>
          <colorvec_file>$colorvec</colorvec_file>
          <soln_files>
EOF
  $self->print_elem_arr($prop_files);

  print OUT <<EOF;
          </soln_files>
          <peram_files>
EOF
  $self->print_elem_arr($peram_files); 

  print OUT <<EOF;
          </peram_files>
          <hadron_op_file>$hadron_node_out_db</hadron_op_file>
          <hadron_node_xmls>
            <elem>$hadron_node_xml</elem>
          </hadron_node_xmls>
        </NamedObject>

      </elem>
EOF
}



sub write_harom_hadron_node_xml
{
  my $self = shift; 
  my $file =$self->get_harom_hadron_node_xml_ini_string();

  if (-e $file)
  { 
    unlink ($file) || die($_); 
  }

  open OUT , ">" , $file;


  $self->print_harom_header(); 
  $self->print_distillution_noise();
  $self->print_harom_hadron_node_dist(); 
  $self->print_harom_trailer(); 

  close OUT;

  return $file;  
}



1;
