#!/usr/bin/perl

use strict; 
use warnings;
use POSIX; 

die("usage: $0 <FF_num> " ) unless $#ARGV == 0;

my @dat = `ls *F_${ARGV[0]}* | grep jack | grep -v fit`;

chomp @dat; 

my @stems = ();

foreach (@dat)
{
  my ($stem,$trash) = split(/\./,$_);
  push @stems , $stem;
}


#print @stems;

my %Q2 = ();


foreach my $stem (@stems)
{
  my ($a,$b,$q2,$c,$d) = split(/_/,$stem);
  $q2 =~ s/p/./g;
  $Q2{$q2} = $stem; 
}


my @qs = sort {$a <=> $b} keys %Q2;


#print @qs; 

mkdir ("tmp"); 


foreach my $q (@qs)
{


  my $rn = int(rand(1000)); 

  open GNU , ">" , "tmp/plot_$rn.gp"; 

  my $stem = $Q2{$q};
  my $cmdr = "calcbc \' real ( ${stem}.jack ) \' > tmp/${stem}.dat.real";
  my $cmdi = "calcbc \' imag ( ${stem}.jack ) \' > tmp/${stem}.dat.imag";

  system ( $cmdr ) == 0 || die($_);
  system ( $cmdi ) == 0 || die($_); 

  my $ffit = `calc ${stem}_fit.jack`;
  my @fit = split /\s+/ , $ffit;

#  print "foobar" . $fit[0] . " " . $fit[1] . " " . $fit[2] . "\n";
  my $high = $fit[1] + $fit[2];
  my $low = $fit[1] - $fit[2];


  print GNU "set label 1 \"$q\" at graph 0.45,0.9 font \',8\' \n";
  print GNU "f_high(x) = $high \n";
  print GNU "f_low(x) = $low \n";
#  print GNU "set style fill transparent solid 0.5 noborder \n";
  print GNU "plot \'tmp/${stem}.dat.imag\' using 1:2:3 w yerr title \'imag\' lc rgb \"blue\", \\\n";
  print GNU " \'+\' using 1:(f_high(\$1)):(f_low(\$1)) with filledcurve title \'fit\' lc rgb \"forest-green\" ,  \\\n";
  print GNU " \'tmp/${stem}.dat.real\' using 1:2:3 w yerr title \'real\' lc rgb \"red\" \n";


  close GNU;


  system ("gnuplot -persist tmp/plot_$rn.gp") ;

}




