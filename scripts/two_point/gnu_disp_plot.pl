#!/usr/bin/perl -w

use strict;



die ("usage: $0 <disp.list> ") unless $#ARGV == 0;

my $file = $ARGV[0]; 

open IN , "<" , $file;
my @lines = <IN>;
close IN;



my %hash = ();

foreach my $line (@lines)
{

  my ($pp,$EE) = split /\s+/ , $line;
#  print "$pp \n $EE \n";

  my $calced = `calc $EE | awk \'{print \$2 " " \$3 }\' `;
#  print $calced . "\n";

  if(exists $hash{$pp} )
  {
    push @{$hash{$pp}} , $calced;
  }
  else
  {
    $hash{$pp} = [$calced];
  }

}


my @keys = sort{$a <=> $b} keys %hash;

mkdir "tmp" unless (-e "tmp");
my $rn = int(rand(1000)); 
my $dat = "tmp/dat_${rn}.dat";
my $plot = "tmp/pt_${rn}.gp";

open OUT , ">" , $dat;
foreach my $key (@keys)
{
  my @vals = @{$hash{$key}};

  foreach my $val (@vals)
  {
    print OUT $key . " " . $val
  }
}

close OUT;

open PLOT , ">" , $plot;
print PLOT "plot \'$dat\' using 1:2:3 w yerr";
close PLOT;

system("gnuplot -p $plot");
