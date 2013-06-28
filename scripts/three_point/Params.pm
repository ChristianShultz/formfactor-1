#!/usr/bin/perl

# make a params package since perl w/o using strict is scary 
package Params;
use strict;

sub new
{
  my $invocant = shift;
  my $class = ref($invocant) || $invocant; 
  my $self =  {
    STEM => "szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per",
    VERSION => 3,
    DIAGNOSTIC => 10,
    SEQNO => undef, 
    NUMVECS => 64,
    NTCORR => 32,
    TSOURCE => 0,
    L_T =>128,
    L_S => 16,
    DECAYDIR => 3,
    CONVERTUDTOL => "true", 
    UMASS => "U-0.0743",
    DMASS => "D-0.0743",
    SMASS => "S-0.0743",    
    LUSTREDIR => "/lustre/volatile/Spectrum/Clover/NF2+1",
    CACHEDIR => "/lustre/cache/Spectrum/Clover/NF2+1",
    WORKDIR => "/work/JLabLQCD/LHPC/Spectrum/Clover/NF2+1",
    HADRONNODEXML => "hadron_node.xml",
    GENGRAPHEXE => "redstar_gen_graph",
    NPTEXE => "redstar_npt",
    MESONEXE => "meson_hadron_node",
    NPOINTLIST => "npt.list.xml",
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


sub nt_corr {
  my $self = shift; 
  if(@_) {$self->{NTCORR} = shift}
  return $self->{NTCORR};
}


sub t_source {
  my $self = shift; 
  if(@_) {$self->{TSOURCE} = shift}
  return $self->{TSOURCE};
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

sub hadron_node_xmls
{
  my $self = shift;
  if(@_) {$self->{HADRONNODEXML} = shift}
  return $self->{HADRONNODEXML};
}

sub gen_graph_exe {
  my $self = shift; 
  if(@_) {$self->{GENGRAPHEXE} = shift}
  return $self->{GENGRAPHEXE};
}


sub npt_exe {
  my $self = shift; 
  if(@_) {$self->{NPTEXE} = shift}
  return $self->{NPTEXE};
}


sub meson_exe {
  my $self = shift; 
  if(@_) {$self->{MESONEXE} = shift}
  return $self->{MESONEXE};
}

sub npt_list{
  my $self = shift;
  if(@_) {$self->{NPOINTLIST} = shift}
  return $self->{NPOINTLIST};
}

#
# More methods.. no args
#



sub ensemble{
  my $self = shift; 
  return $self->{STEM};
}


sub prop_dbs{
  my $self = shift; 
  my @prop_dbs; 
  my $cache_dir = $self->cache_dir();
  my $stem = $self->stem();
  my $seqno = $self->seqno();
  push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_0-124_inc4.sdb${seqno}");
  push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_1-125_inc4.sdb${seqno}");
  push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_2-126_inc4.sdb${seqno}");
  push(@prop_dbs, "$cache_dir/${stem}/prop_db/${stem}.prop.t0_3-127_inc4.sdb${seqno}");

  return \@prop_dbs; 
}

sub meson_dbs
{
  my $self = shift;
  my @mdbs;
  my $stem = $self->cache_dir() . "/" . $self->stem() ."/meson_db_dispmom/" . $self->stem() ;
  push @mdbs , "${stem}.meson.colorvec.sdb" . $self->seqno(); 
  return \@mdbs; 
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

sub hadron_node_graph
{
  my $self = shift; 
  my $db = "graph.sdb" . $self->seqno();
  return $db; 
}

sub hadron_node_output_db
{
  my $self = shift; 
  my $db = "hadron_node.sdb".$self->seqno();
  return $db; 
}

sub output_db
{
  my $self = shift; 
  my $db = "output.npt.sdb" . $self->seqno();
  return $db; 
}

sub print_meson_hadron_node_Param
{
  my $self = shift; 

  print OUT <<EOF;
  <Param>
    <version>$self->{VERSION}</version>
    <num_vecs>$self->{NUMVECS}</num_vecs>
    <Nt_corr>$self->{NTCORR}</Nt_corr>
    <t_origin>$self->{TSOURCE}</t_origin>

    <FlavorToMass>
      <elem>
        <flavor>l</flavor>
        <mass>$self->{UMASS}</mass>
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
      <elem>$self->{HADRONNODEXML}</elem>
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

sub write_meson_hadron_node_xml
{
  my $self = shift; 
  my $file = "meson_hadron_node.xml" . $self->seqno(); 

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


sub print_redstar_npoint_list
{
  my $self = shift; 
  my $file = $self->npt_list(); 
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
  my $had_graph = $self->hadron_node_graph(); 
  my $had_node_out_db = $self->hadron_node_output_db();
  my $out_db = $self->output_db();  
  print OUT <<EOF;
  <DBFiles>
    <proj_op_xmls></proj_op_xmls>
    <hadron_node_xml>$self->{HADRONNODEXML}</hadron_node_xml>
    <hadron_npt_graph_db>$had_graph</hadron_npt_grap_db>
    <hadron_node_dbs><elem>$had_node_out_db</elem></hadron_node_dbs>
    <output_db>$out_db</output_db>
  </DBFiles>
EOF
}

sub write_redstar_xml
{
  my $self = shift;
  my $file = "redstar_xmlini.xml" . $self->seqno();

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


1;
