#!/usr/bin/perl -w

use strict;



die ("usage: $0 <f.list> <mass> <aniso> <L_s> <output=optional>") unless $#ARGV >= 3;

my $file = $ARGV[0]; 
my $mass = $ARGV[1];
my $xi = $ARGV[2]; 
my $L_s = $ARGV[3]; 
my $outgraph = "output.eps";

if ( $#ARGV == 4 ) 
{
  $outgraph = $ARGV[4]; 
}

open IN , "<" , $file;
my @lines = <IN>;
close IN;


my $mass2 = $mass*$mass; 
my $unit = 2.*3.14159/($xi *$L_s); 

my %hash = ();

foreach my $line (@lines)
{

  my ($pp,$EE) = split /\s+/ , $line;
#  print "$pp \n $EE \n";

  next unless $pp != 0;

  my $denom = $unit*$unit*$pp; 


  my $calced = `calcbc ' ( $EE * $EE - $mass2 ) / ( $denom )' | awk \'{print \$2 " " \$3 }\' `;

  
  
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
print PLOT "g(x) = c*c \n"; 

print PLOT "fit g(x) \'$dat\' using 1:2 via c \n ";
print PLOT "plot \'$dat\' using 1:2:3 w yerr ti '(E^2(p^2) - m^2 )/p^2', \\\n";
print PLOT "g(x) ti sprintf(\"g-> c=%.4f\", c)  \n";
print PLOT "set xr [GPVAL_DATA_X_MIN - 0.1 : GPVAL_DATA_X_MAX + 0.1] \n";
print PLOT "replot \n";
print PLOT "set terminal postscript eps color \n";
print PLOT "set out \'$outgraph\' \n";
print PLOT "replot \n";
print PLOT "set term pop \n";
close PLOT;

system("gnuplot -p $plot");
