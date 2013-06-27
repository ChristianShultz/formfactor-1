#!/usr/bin/perl

die "usage: $0 <weights.particle.list>" unless $#ARGV == 0; 

use strict;

# some global switch to turn off/on unit normalization
# NB: only undef is false since perl is dumb
our $donorm = undef;  
our $saveplot = "true";

# the big list file
my $file = $ARGV[0];

# read data
open IN , "<" , $file || die $_; 
my @lines = <IN>;
close IN; 
chomp @lines; 

# init a hash
my %hash = ();


# process the data.. key is some optimized op name, dat is the interpolating field and weight
foreach my $line (@lines)
{
  my ($proj_op, $op, $weight) = split(/ /, $line);
  my $key = $proj_op;
  my $dat = $op . " " . $weight;
  push @{ $hash{$key} } , $dat if exists $hash{$key};  
  $hash{$key} = [$dat] unless exists $hash{$key};
}

# our standard rng naming and a tmp dir
my $tmp = "tmp";
my $rn = int(rand(1000));
mkdir $tmp;

# get the set of optimized operators
my @keys = keys %hash;

# make it do something pretty
print "extracting.. \n"; 

# the set of plots (one for each optimized operator) 
my @plots = (); 

# loop over each optimized operator and generate the plot 
foreach my $key (@keys)
{
  print $key . "\n";
  push @plots, make_gp($rn,$key,$hash{$key});  
}

# loop and run gnuplot
foreach my $plot (@plots)
{
  system("gnuplot -persist $plot"); 
}


# make the plot script
#######################################################
#######################################################
#######################################################
sub make_gp
{
  my ($rn,$tag, $aref) = @_; 

  my @dat = @{$aref};

  my ($opr,$sgnr) = normalize(\@dat); 

my %ops = %{$opr};
my %sgns = %{$sgnr};
my @keys = keys %ops; 

# a sort based on spins
my @j0 = ();
my @j1 = ();
my @j2 = ();
my @j3 = ();
my @j4 = (); 

# assuming standard operator names __J#
foreach my $key (@keys)
{
  push @j0, $key if $key =~ m/__J0/;
  push @j1, $key if $key =~ m/__J1/;
  push @j2, $key if $key =~ m/__J2/;
  push @j3, $key if $key =~ m/__J3/;
  push @j4, $key if $key =~ m/__J4/;
}

# set up a tag for this running of the perl script
my $file = "tmp/".$tag . "_" .$rn;


# dump the data into a file for gnuplot
my $fout = $file . "_" . ".dat";
open DATOUT , ">" , $fout;

my $xp = 0.; 
my $j0f = dump_data($file,0,\@j0,\%ops,\%sgns);
my $j1f = dump_data($file,1,\@j1,\%ops,\%sgns); 
my $j2f = dump_data($file,2,\@j2,\%ops,\%sgns); 
my $j3f = dump_data($file,3,\@j3,\%ops,\%sgns); 
my $j4f = dump_data($file,4,\@j4,\%ops,\%sgns); 

close DATOUT;


# make the plot
my $plot = $file . ".gp";

open OUT , ">" , $plot; 

print OUT "set xtics nomirror rotate by -45 font \",8\" \n";
print OUT "set key outside \n"; 
print OUT "set style data histogram \n";
print OUT "set style fill solid border -1 \n";
print OUT "set title \"$tag\" \n";
print OUT "unset key \n";
print OUT "set rmargin at screen 0.85 \n";
print OUT "set yrange [0:1.05] \n";
print OUT "set object 1 rect from screen 0.90,0.85 to screen 0.95,0.90 fc lt 0 \n";
print OUT "set object 2 rect from screen 0.90,0.75 to screen 0.95,0.80 fc lt 1 \n";
print OUT "set object 3 rect from screen 0.90,0.65 to screen 0.95,0.70 fc lt 2 \n";
print OUT "set object 4 rect from screen 0.90,0.55 to screen 0.95,0.60 fc lt 3 \n";
print OUT "set object 5 rect from screen 0.90,0.45 to screen 0.95,0.50 fc lt 4 \n";
print OUT "set label 1 at screen 0.915,0.87 \"J0\" front \n";
print OUT "set label 2 at screen 0.915,0.77 \"J1\" front \n";
print OUT "set label 3 at screen 0.915,0.67 \"J2\" front \n";
print OUT "set label 4 at screen 0.915,0.57 \"J3\" front \n";
print OUT "set label 5 at screen 0.915,0.47 \"J4\" front \n";
print OUT "plot \"$fout\" u (column(0)):2:(0.5):(\$3):xtic(1) w boxes lc variable, \\\n" ;
#                                #xval:ydata:boxwidth:color_index:xtic_labels
print OUT "\"$fout\" u (column(0)):(\$2 + 0.02):4 with labels";


if($saveplot)
{
# looks like crap
#  my $png = $file . ".png";
# print OUT "\nset output \"$png\" \n";
#  print OUT "set term png\n";

  my $ps = $tag . ".eps";
  print OUT "\nset terminal postscript eps \n";
  print OUT "set output \"$ps\" \n";
  print OUT "replot"

}


#print OUT "plot \"$fout\" using xticlabels(1):2:3   \n"; 

close OUT; 


return $plot; 
}


# unit normalize and keep track of the sign
#######################################################
#######################################################
#######################################################
sub normalize
{
  my ($ref) = @_; 
  my @in = @{$ref};
  my @out = (); 
  my $largest = 0; 
  my %ops = (); 
  my %sign = (); 
  foreach my $elem (@in)
  {
    my ($op, $val) = split(/ /, $elem); 

    if($val < 0)
    {
      $ops{$op} = -$val;
      $sign{$op} = "-";
    }
    else
    {
      $ops{$op} = $val;
      $sign{$op} = "+";
    }

    if ($ops{$op} > $largest)
    {
      $largest = $ops{$op};
    }

  }

  my @keys = keys %ops;

  if($donorm)
  {
    foreach my $key (@keys)
    {
      my $foo = $ops{$key};
      $ops{$key} = $foo/$largest; 
    }
  }

  return (\%ops, \%sign); 
}

# dump the data into something gnuplot can read
#######################################################
#######################################################
#######################################################
sub dump_data
{
  my ($file,$j,$ky,$opv,$sg) = @_;

  my @keys = @{$ky};
  my %ops = %{$opv};
  my %sgn = %{$sg};

  foreach my $key (@keys)
  {
    print DATOUT $key . " " . $ops{$key} . " " . $j  . " " . $sgn{$key} . "\n"; 
  }

  return $#keys; 
}


