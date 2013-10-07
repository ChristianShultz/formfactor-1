#
# 
# WORK SCRIPT NEEDED BY OPperams.pl
#
#

use strict; 
use File::Basename;
use File::Copy;

sub rep_list
{
  my ($spin , $parity) = @_; 

  print "generating irrep list for spin $spin , parity $parity \n";

  my %fcns = (); 

  $fcns{0} = \&rep_list_spin_0;
  $fcns{1} = \&rep_list_spin_1;

  my $ref = "";
  if ( exists $fcns{$spin} )
  {
    $ref = $fcns{$spin}; 
  }
  else
  {
    my @kys = keys %fcns; 
    die( "error, supported spins are " . @kys . "\n" );
  }

  return &$ref($parity);
}



sub rep_list_spin_0
{
  my $parity = $_[0]; 

  my %ret = (); 
  my $rep = "";

  # etap = parity
  if ( $parity eq "p" )
  {
    $rep = "A1";
  }
  elsif ( $parity eq "m" )
  {
    $rep = "A2";
  }
  else
  {
    die ( " parity must be either p or m , parity : $parity ");
  }

  my %rest = (); 
  $rest{ "A1" } = "";

  my %D2 = (); 
  $D2{ "D2" . $rep } = "H0";

  my %D3 = (); 
  $D3{ "D3" . $rep } = "H0";

  my %D4 = (); 
  $D4{ "D4" . $rep } = "H0";


  $ret { "p000" } = \%rest; 
  $ret { "p100" } = \%D4; 
  $ret { "p110" } = \%D2; 
  $ret { "p111" } = \%D3; 
  $ret { "p200" } = \%D4; 

  return \%ret; 

}



sub rep_list_spin_1
{
  my $parity = $_[0]; 

  my %ret = (); 
  my $rep = "";

  # etap = -parity
  if ( $parity eq "p" )
  {
    $rep = "A2";
  }
  elsif ( $parity eq "m" )
  {
    $rep = "A1";
  }
  else
  {
    die ( " parity must be either p or m , parity : $parity ");
  }

  my %rest = (); 
  $rest{ "T1" } = "";

  my %D2 = (); 
  $D2{ "D2" . $rep } = "H0";
  $D2{ "D2B1" } = "H1";
  $D2{ "D2B2" } = "H1";

  my %D3 = (); 
  $D3{ "D3" . $rep } = "H0";
  $D3{ "D3E2" } = "H1";

  my %D4 = (); 
  $D4{ "D4" . $rep } = "H0";
  $D4{ "D4E2" } = "H1";


  $ret { "p000" } = \%rest; 
  $ret { "p100" } = \%D4; 
  $ret { "p110" } = \%D2; 
  $ret { "p111" } = \%D3; 
  $ret { "p200" } = \%D4; 

  return \%ret; 

}




sub print_header_perl
{

  print OUT <<EOF
#!/usr/bin/perl -w
#
# this is some auto generated perl skeleton
# and should not be assumed to be correct or
# work in the way you think it should
#


use strict; 
use OPparams; 
use File::Basename;

my \$basedir = dirname(\$0); 

require "\${basedir}/OPparams.work.pl";

my \@all_ops = (); 

EOF

}



sub print_single_op
{
  my ($pid, $rep , $rep_stem , $mom , $twoI_z, $ncfg , $ensemble , $recon_dir ) = @_;  
  my $p = "\$${pid}_${mom}_${rep_stem}";


  print OUT <<EOF;

  my $p = OPparams->new();
  $p->pid(\"$pid\");
  $p->irrep(\"$rep\");
  $p->irrep_stem(\"$rep_stem\");
  $p->mom(\"$mom\"); 
  $p->twoI_z($twoI_z); 
  $p->ncfg($ncfg);
  $p->ensemble(\"$ensemble\"); 
  $p->recon_dir(\"$recon_dir\");
  $p->t0(-1);
  $p->state(-1); 
  $p->tz(-1);
  $p->phaser(1.);
EOF

  return $p; 
}


sub print_ops
{

  my ($pid, $twoI_z, $ncfg , $ensemble , $recon_dir , $hashr ) = @_;  


  my %hash = %{ $hashr };
  my @kys = keys %hash; 

  @kys = sort { $a cmp $b } @kys; 

  my @all_ops = (); 

  foreach my $k (@kys)
  {
    #    print "working on $k \n";

    my $rephashr = $hash{$k};
    my %rephash = %{$rephashr};

    my @repkys = keys %rephash;

    @repkys = sort { $a cmp $b } @repkys;

    foreach my $rk (@repkys)
    {
      print "mom = $k rep = $rk \n"; 
      push @all_ops , &print_single_op($pid , $rk , $rephash{$rk}.$rk , $k, $twoI_z , $ncfg , $ensemble , $recon_dir ); 
    }
  } 


  print OUT <<EOF;


  my \@all_operators = (); 
EOF



  foreach my $p (@all_ops)
  {
    print OUT <<EOF;
  push \@all_operators , $p; 
EOF
  }

}


sub print_trailer_perl
{

  my $pid = shift; 

  print OUT <<EOF;


  my \@extracts = (); 
  foreach my \$op (\@all_operators)
  {
    push \@extracts,  &run_extract_all_v_coeffs_svd(\$op); 
  }

  my \$listfile = &convert_proj_to_xml(\\\@extracts,\"$pid\"); 

  &make_proj_plots(\$listfile);

  &write_radmat_xml(\\\@all_operators); 
EOF


}


sub write_radmat_xml
{
  my $ref = shift; 
  my @all_ops = @{$ref};
  my %hash = (); 
  my $dir = `pwd`;
  chomp $dir;  
  my $rdir = $dir . "/radmat";
  mkdir $rdir unless -d $rdir;

  #move down
  chdir $rdir; 

  # a month later I am sorry for writing this..eww
  foreach my $op (@all_ops)
  {
    my %h = %{ $op->write_mass_overlap_xml() };

    # hash slicing -- keys are op names which are unique
    # so this is a safe but cryptic operation 
    @hash{ keys %h } = values %h ; 
  }

  # make some lists so we can do that radmat thing later 
  open OUT , ">" , "mass_overlap_xml_files.list";
  foreach my $k (keys %hash)
  {
    print OUT $k . "\n";
  }
  close OUT; 

  open OUT , ">" , "mass_overlap_sdbs.list"; 
  foreach my $v (values %hash) 
  {
    print OUT $v . "\n";
  }
  close OUT; 

  #move back out
  chdir $dir; 
}



sub copy_file
{
  my ($f,$dest) = @_; 

  die ( "ERROR : $f doesnot exist \n" ) unless ( -f $f );


  if ( -f $dest ) 
  {
    unlink $dest; 
  }

  copy($f,$dest) || die ("could not copy $f to $dest \n $_");  

}



sub run_extract_all_v_coeffs_svd
{
  my $op = shift;
  my $base = `pwd`;
  chomp $base;  
  my $destdir = $op->recon_dir(); 


  my $loc = "/u/home/shultz/optimized_operators/";
  my $exe = $loc ."extract_all_v_coeffs_svd.pl";
  my $t0 = $op->t0(); 
  my $tz = $op->tz(); 
  my $state = $op->state(); 
  my $opslistfile = "ops_phases";
  my $opname = $op->op_name(); 

  my $outfile = $opname . ".list"; 

  my $run = "${exe} $t0 $tz $state $opslistfile $opname > $outfile";

  chdir $destdir || die ( $_ ); 

  my $script =  "${opname}.extract_command.csh"; 
  open OP , ">" , $script; 
  print OP "#/bin/tcsh \n";
  print OP $run ." \n"; 
  close OP;

  # for some reason they go 0644 for me..
  chmod 0755 , $exe;
  chmod 0755 , $script; 

  system ( " ./${script} " ) == 0 || die ( "some problem for $run \n in $destdir \n" ) ; 


  chdir $base || die ( $_ );


  copy_file($destdir."/".$outfile , $outfile);

  return $outfile; 
}



sub convert_proj_to_xml
{
  my ( $ref , $id ) = @_; 
  my @lists = @{$ref};

  my $listf = "weights." . $id . ".list";
  my $xmlf = "weights." . $id . ".xml";

  unlink $listf unless ! -f $listf; 
  unlink $xmlf unless ! -f $xmlf; 

  open OUT , ">" , $listf; 

  foreach my $l (@lists)
  {
    open IN , "<" , $l; 
    print OUT while <IN>; 
    close IN; 
  } 

  close OUT; 

## THIS IS OLD
#  my $exe = "convert_proj_list_to_xml.pl"; 
#  die ( "$exe not present" ) unless -f $exe; 
#  
#  system (" ./${exe} $listf $xmlf " ) == 0 || die ($_) ;

  my $loc = "/u/home/shultz/optimized_operators/";
  my $exe = $loc . "convert_proj_list_to_irrep_op_xml.pl";
  die ("$exe not present") unless -f $exe;
  system (" cat $listf | ${exe} > $xmlf ") == 0 || die($_); 

  return $listf;  
}

sub make_proj_plots
{
  my $file = shift; 
  my $loc = "/u/home/shultz/optimized_operators/";
  my $exe = $loc . "plot_proj_op_coeffs.pl"; 

  die ( "$exe not present" ) unless -f $exe; 
  die ( "input not present : $file " ) unless -f $file; 

  system ( " ${exe} $file " ) == 0 || die ($_); 
}

1;  # for require 
