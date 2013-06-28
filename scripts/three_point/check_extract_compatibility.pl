#!/usr/bin/perl -w


use strict;

my @keys = ("p000","p100","p110","p111","p200");

my %vals = ();


# get the irrep list and energy
# vals key is momentum , value is a hash where irrep->energy

foreach my $mom (@keys)
{
  my $num_elem = `ls *$mom*$mom* | grep _E | wc -l`;

  if ($num_elem > 1)
  {
    my $elems = `ls *$mom*$mom* | grep _E | xargs`;

    my @irreps = split(/\s+/,$elems);

    my %rep_E = ();
    foreach my $rep (@irreps)
    {
      $rep_E{$rep} = `calc $rep | awk '{print \$2}'`;
    }
    $vals{$mom} = \%rep_E;
  }
}



my @kys = keys %vals; 
foreach my $key (@kys)
{
  my %hash = %{ $vals{$key} };

  my @reps = keys %hash; 
# loop every pair looking for inconsistency   
  for my $i (0..$#reps)
  {
    for my $j ($i+1..$#reps)
    {
      my $ki = $reps[$i];
      my $kj = $reps[$j];
      
        print $ki . " vs " . $kj . " percent diff = " . percent_diff($hash{$ki},$hash{$kj}) . "\n" ; 

      print $ki . " is inconsistent with " . $kj . " percent diff = " . percent_diff($hash{$ki},$hash{$kj}) unless percent_diff($hash{$ki},$hash{$kj}) < 0.05;
    }
  }

}


sub percent_diff
{
  my ($E1,$E2) = @_; 

  my $E = $E1 > $E2 ? $E1 : $E2; 


  return abs($E1 - $E2) / $E; 
}
