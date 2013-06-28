#!/usr/bin/perl

use strict; 
use warnings;
use POSIX; 

# usage 
die("usage: $0 <file.jack> <make_ratio present = true>") unless $#ARGV == 0; 


# params
my $file = $ARGV[0];
my $dorat = $ARGV[1];

print $file;


# trash dir
mkdir ("tmp"); 

my $rn = int(rand(1000));

mkdir ("tmp");


#### DIEGO -- the next set of lines assume that you have an ensemble and not a single configuration 
#             you need to do some sort of hack to get the code to shove only the real/imag parts 
#             into the files 


# make data command
my $cmdr = "calcbc \' real ( $file ) \' > tmp/${file}.real";
my $cmdi = "calcbc \' imag ( $file ) \' > tmp/${file}.imag";
my $cmdrat = "calcbc \' real( $file) / imag ( $file ) \' > tmp/${file}.ratio";

# execute make data command
system ( $cmdr ) == 0 || die($_);
system ( $cmdi ) == 0 || die($_);
system ( $cmdrat ) == 0 || die($_);


# shove this junk into a gnuplot file
open GNU , ">" , "tmp/plot_${rn}.gp"; 
#  print GNU "set label 1 \"$q\" at graph 0.35,0.9 font \',8\' \n";
print GNU "set title \"";


# split up the name to be 80 char then new line then next 80 char...
my @chars = map substr($file, $_, 1), 0 .. length($file) -1; 
my $idx = 0;
foreach my $char (@chars) 
{ 
  if($idx == 80)
  {
    print GNU "\\n";
    $idx = 0; 
  }
  print GNU $char;
  $idx++; 
}

print GNU "\"\n";
print GNU "plot \'tmp/${file}.real\' using 1:2:3  w yerr title \'real\', \\\n";
print GNU " \'tmp/${file}.imag\' using (\$1 + 0.5):2:3 lc 3 w yerr title \'imag\' \n";
print GNU "unset multiplot";
close GNU;


# run the gnuplot
system ("gnuplot -persist tmp/plot_${rn}.gp") ;



# some crap that you don't need 

if ($dorat)
{
  open GNUR , ">" , "tmp/plot_${rn}.rat.gp";

  print GNUR "set title \"";

  $idx = 0;
  foreach my $char (@chars) 
  { 
    if($idx == 80)
    {
      print GNUR "\\n";
      $idx = 0; 
    }
    print GNUR $char;
    $idx++; 
  }

  print GNUR "\"\n";
  print GNUR "plot \'tmp/${file}.ratio\' using 1:2 title \'ratio \' \n";
  close GNUR;
  system ("gnuplot -persist tmp/plot_${rn}.rat.gp");

}


#system ("rm tmp/plot_$rn.gp tmp/${file}.real tmp/${file}.imag tmp/plot_${rn}.rat.gp "); 


