#!/usr/bin/perl

use strict; 
use warnings;
use POSIX; 



my @dat = `ls | grep jack`;

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

  my $sstem = substr $stem , 14;
  $sstem =~ s/_//g;
#  print $sstem;

  $Q2{$sstem} = $stem; 
}


my @qs = sort {$a cmp $b} keys %Q2;

my $num = ceil(sqrt($#qs)); 


#print @qs; 

mkdir ("tmp"); 

my $rn = int(rand(1000));


if ($num < 5)
{

  open GNU , ">" , "tmp/plot_$rn.gp"; 

  print GNU "set multiplot layout $num , $num rowsfirst \n ";

  foreach my $q (@qs)
  {
    my $stem = $Q2{$q};
    my $cmdr = "calcbc \' real ( ${stem}.jack ) \' > tmp/${stem}.dat.real";
    my $cmdi = "calcbc \' imag ( ${stem}.jack ) \' > tmp/${stem}.dat.imag";

    system ( $cmdr ) == 0 || die($_);
    system ( $cmdi ) == 0 || die($_); 


#  print GNU "set label 1 \"$q\" at graph 0.35,0.9 font \',8\' \n";
    print GNU "set title \"$q\" \n";
    print GNU "plot \'tmp/${stem}.dat.real\' using 1:2:3  w yerr title \'real\', \\\n";
    print GNU " \'tmp/${stem}.dat.imag\' using 1:2:3 lc 3 w yerr title \'imag\' \n";
  }


  print GNU "unset multiplot";
  close GNU;


  system ("gnuplot -persist tmp/plot_$rn.gp") ;
}


if ($num > 5)
{

  my $chunksz = 16;
  $num = sqrt($chunksz);
  my $nchunk = ceil( $#qs / $chunksz ); 

  for my $chunk (0..$nchunk)
  {

    open GNU , ">" , "tmp/plot_$rn.gp"; 
    print GNU "set multiplot layout $num , $num rowsfirst \n ";


    for my $work (0..$chunksz)
    {
      my $q = $qs[$chunk * $chunksz + $work];
      my $stem = $Q2{$q};
      my $cmdr = "calcbc \' real ( ${stem}.jack ) \' > tmp/${stem}.dat.real";
      my $cmdi = "calcbc \' imag ( ${stem}.jack ) \' > tmp/${stem}.dat.imag";

      system ( $cmdr ) == 0 || die($_);
      system ( $cmdi ) == 0 || die($_); 


#  print GNU "set label 1 \"$q\" at graph 0.35,0.9 font \',8\' \n";
      print GNU "set title \"$q\" \n";
      print GNU "plot \'tmp/${stem}.dat.real\' using 1:2:3  w yerr title \'real\', \\\n";
      print GNU " \'tmp/${stem}.dat.imag\' using 1:2:3 lc 3 w yerr title \'imag\' \n";


    }

    print GNU "unset multiplot";
    close GNU;


    system ("gnuplot -persist tmp/plot_$rn.gp") ;


  }

}


