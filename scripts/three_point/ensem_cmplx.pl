#!/usr/bin/perl 

use strict;
use warnings; 



die "$0 : usage <real> <cmplx> <bins> <time> <type 0=real, 1=complex> <space> <#cols>" unless $#ARGV == 6;

my $file = "ensem_num.jack"; 

open OUT , ">" , $file || die($_);


print OUT "$ARGV[2] $ARGV[3] $ARGV[4] $ARGV[5] $ARGV[6] \n";


for my $bn (0..$ARGV[2])
{
  for my $time (0..$ARGV[3] -1)
  {
    print OUT $time . " " . $ARGV[0];
    print OUT " " . $ARGV[1] unless ($ARGV[4] == "0");
    print OUT "\n";
  }
}

close OUT;




