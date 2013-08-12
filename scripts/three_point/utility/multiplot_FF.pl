#!/usr/bin/perl

use strict; 
use warnings;
use POSIX; 

die("usage: $0 <FF_num> " ) unless $#ARGV == 0;

my @dat = `ls *F_${ARGV[0]}* | grep jack`;

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

my $num = ceil(sqrt($#qs)); 


#print @qs; 

mkdir ("tmp"); 

my $rn = int(rand(1000)); 

open GNU , ">" , "tmp/plot_$rn.gp"; 

print GNU "set multiplot layout $num , $num rowsfirst \n ";

foreach my $q (@qs)
{
  my $stem = $Q2{$q};
  my $cmdr = "calcbc \' real ( ${stem}.jack ) \' > tmp/${stem}.dat.real";
  my $cmdi = "calcbc \' imag ( ${stem}.jack ) \' > tmp/${stem}.dat.imag";

  system ( $cmdr ) == 0 || die($_);
  system ( $cmdi ) == 0 || die($_); 


  print GNU "set label 1 \"$q\" at graph 0.45,0.9 font \',8\' \n";
  print GNU "plot \'tmp/${stem}.dat.real\' using 1:2:3 w yerr title \'real\', \\\n";
  print GNU " \'tmp/${stem}.dat.imag\' using 1:2:3 w yerr title \'imag\' \n";
}


print GNU "unset multiplot";
close GNU;


system ("gnuplot -persist tmp/plot_$rn.gp") ;

