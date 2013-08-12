#!/usr/bin/perl -w


BEGIN{
  push (@INC, "/Users/shultz/jlab/scripts");
}

use strict; 
use findPhase; 

my @qs = split (/\s+/,  ` ls -al | grep ^d | awk '{ print \$9 }' | grep Q2_ | xargs ` );

my @lines;
my $gsum = 0; 
my $zsum = 0; 

foreach my $q2 (@qs)
{
  chomp $q2;

  my $q2dir = $q2 . "/llsq/" ;

  my $gnum = -1; 

  if ( ! -d $q2dir )
  {
    print " warning $q2 has no llsq dir.. outside analytic region? \n";
    $gnum = 0; 
  }

  $gnum = `ls $q2dir/*.jack 2>/dev/null | wc -l | awk '{ print \$1}'` unless $gnum == 0;

  my $f = "foobar";
  $f = `head -1 $q2dir/row_index_to_continuum_elem.txt | awk '{ print \$2 }'` unless $gnum == 0;

  my $pf = "undef"; 
  my $pi = "undef";

  if( $gnum != 0)
  {
    my $ph = findPhase->new();
    $ph->id($f); 
    $ph->parse_id(); 
    $pf = $ph->psnk_str(); 
    $pi = $ph->psrc_str(); 
  }

  my $zero_dir = $q2dir . "zeroed_matrix_elems/";

  my $znum = -1; 

  if ( ! -d $zero_dir )
  {
    print " warning $q2 has no zeroed llsq dir.. outside analytic region? \n";
    $znum = 0; 
  }


  $znum =  `ls $zero_dir/*.jack 2>/dev/null | wc -l  | awk '{print \$1}'` unless $znum == 0;


  chomp $gnum; 
  chomp $znum;

  $gsum = $gsum + $gnum; 
  $zsum = $zsum + $znum;

  my ($trash,$qq) = split /_/ , $q2;
  my @line = ($qq, $pf, $pi, $gnum , $znum) ;
  push @lines , \@line;
}


my %hash = (); 

foreach my $line (@lines) 
{
my @l = @{$line};
$hash{$l[0]} = \@l; 
}

my @kys = sort {$a <=> $b} keys %hash; 

open OUT , ">" , "table.tex.txt";

print OUT <<EOF;
\\begin{tabular}{|c|c|c|c|c|}
\\hline
Q^2 & p_f & p_i & nele & nzero \\\\ \\hline
EOF

foreach my $ky (@kys)
{
  my @l = @{ $hash{$ky} };

  my $ee = pop @l; 
  foreach my $e (@l)
  {
    print OUT $e . " & ";
  }
  print OUT $ee . "\\\\ \\hline \n"; 
}

print OUT "TOTAL & & & $gsum & $zsum \\\\ \\hline \n";

print OUT <<EOF;
\\end{tabular}
EOF

close OUT; 

