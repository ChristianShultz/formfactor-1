#!/usr/bin/perl -w

use strict;



die ("usage: $0 <disp.list> <L_s> <mass=optional> <outfile=optional>") unless $#ARGV >= 1;

my $file = $ARGV[0]; 
my $L_s = $ARGV[1]; 
my $mass = undef; 
my $outgraph = "plot.eps";

if ( $#ARGV >= 2) 
{
  $mass = $ARGV[2]; 
}

if ( $#ARGV == 3 ) 
{
  $outgraph = $ARGV[3];
}

open IN , "<" , $file;
my @lines = <IN>;
close IN;



my %hash = ();

foreach my $line (@lines)
{

  my ($pp,$EE) = split /\s+/ , $line;
#  print "$pp \n $EE \n";

  my $calced = `calcbc ' $EE * $EE ' | awk \'{print \$2 " " \$3 }\' `;
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
print PLOT "t = 6.2831853 \n";
print PLOT "l = $L_s \n";
print PLOT "f = (t/l)*(t/l) \n";

print PLOT "m = $mass\n" unless ! $mass; 
print PLOT "g(x) = m*m + (f/(a*a))*x\n"; 

print PLOT "fit g(x) \'$dat\' using 1:2 via m,a \n" unless $mass;
print PLOT "fit g(x) \'$dat\' using 1:2 via a \n" unless ! $mass;

print PLOT "plot \'$dat\' using 1:2:3 w yerr ti 'E^2(p^2)', \\\n";
print PLOT "g(x) ti sprintf(\"g=m^2 + x*(2pi/Ls/xi)^2, m=%.4f, xi=%.4f\", m, a)  \n";

print PLOT "set terminal postscript eps color \n";
print PLOT "set out \'$outgraph\' \n";
print PLOT "replot \n";
print PLOT "set term pop \n";
close PLOT;

system("gnuplot -p $plot");
