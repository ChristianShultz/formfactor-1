#!/usr/bin/perl

use strict; 
use warnings;
use POSIX; 


die("usage: $0 <mode 1->multiplot, 0->single graphs >") unless ($#ARGV == 0);
my $mode = $ARGV[0];

my @dat = `ls *.jack`;

chomp @dat; 

my @stems = ();

foreach (@dat)
{
  my ($stem,$trash) = split(/\./,$_);
  push @stems , $stem;
}

@stems = sort{do_split_last($a) <=> do_split_last($b)}  @stems;

#print @stems;

my %Q2 = ();

my $num = ceil(sqrt($#dat + 1)); 
my $maxnum = 4;

$num = $maxnum unless ($num < $maxnum);

my @chunks = @{ chunk($num,@stems) };

print "this is stupid $#chunks \n ";

#print @qs; 

mkdir ("tmp"); 

if($mode != 0)
{

  foreach my $this_chunk (@chunks)
  {
    my @these_stems = @{$this_chunk};

    my $rn = int(rand(1000)); 

    open GNU , ">" , "tmp/plot_$rn.gp"; 


    print GNU "set multiplot layout $num , $num rowsfirst \n ";

    foreach my $stem (@these_stems)
    {

      my $cmdr = "calcbc \' real ( ${stem}.jack ) \' > tmp/${stem}.dat.real";
      my $cmdi = "calcbc \' imag ( ${stem}.jack ) \' > tmp/${stem}.dat.imag";

      system ( $cmdr ) == 0 || die($_);
      system ( $cmdi ) == 0 || die($_); 


      print GNU "set label 1 \"$stem\" at graph 0.45,0.9 font \',8\' \n";
      print GNU "plot \'tmp/${stem}.dat.real\' using 1:2:3 w yerr title \'real\', \\\n";
      print GNU " \'tmp/${stem}.dat.imag\' using 1:2:3 w yerr title \'imag\' \n";
    }


    print GNU "unset multiplot";
    close GNU;

    system ("gnuplot -persist tmp/plot_$rn.gp") ;
  }
} 
else
{
  foreach my $stem (@stems)
  {
    my $rn = int(rand(1000)); 

    open GNU , ">" , "tmp/plot_$rn.gp"; 

    my $cmdr = "calcbc \' real ( ${stem}.jack ) \' > tmp/${stem}.dat.real";
    my $cmdi = "calcbc \' imag ( ${stem}.jack ) \' > tmp/${stem}.dat.imag";

    system ( $cmdr ) == 0 || die($_);
    system ( $cmdi ) == 0 || die($_); 


    print GNU "set label 1 \"$stem\" at graph 0.45,0.9 font \',8\' \n";
    print GNU "plot \'tmp/${stem}.dat.real\' using 1:2:3 w yerr title \'real\', \\\n";
    print GNU " \'tmp/${stem}.dat.imag\' using 1:2:3 w yerr title \'imag\' \n";
    close GNU;

    system ("gnuplot -persist tmp/plot_$rn.gp") ;

  }

}


sub chunk
{
  my ($num,@stems) = @_; 
  my $n2 = $num*$num;
  my $nchunk = ceil($#stems/$n2); 
  my @chunks = ();

  for my $i(0..$nchunk-1)
  {
    my @this_chunk = ();

    for my $idx(0..$n2-1)
    {
      push @this_chunk , $stems[$i*$n2 + $idx] unless $i*$n2 + $idx > $#stems;
    }

#    print join(", ",@this_chunk) . "\n";

    push @chunks , \@this_chunk;
  }


  print join(",   ", @chunks);

  return \@chunks; 
}


sub do_split_last
{
  my $foo = $_[0];
  my @bar = split(/_/,$foo);
  return $bar[$#bar];
}
