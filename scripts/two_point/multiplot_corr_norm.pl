#!/usr/bin/perl

use strict;
use warnings; 
use POSIX;


die("usage: $0 <chunksize>") unless ($#ARGV == 0);
my $chunksize = $ARGV[0];


my @dat = `ls *dat*`;
chomp @dat; 


my %stems = ();

foreach (@dat)
{
  my ($stem,$junk) = split(/\.dat_/,$_ ); 

  $stems{$stem} = 0 unless exists $stems{$stem}; 

  $stems{$stem}++; 
}


my @kys = keys %stems; 

foreach (@kys) 
{
  die ("couldn't find everything to make graph for $_") unless $stems{$_} == 3; 
}


mkdir ("tmp"); 

my $num = $#kys + 1; 
my $nchunk = $num/$chunksize; 

for my $chunk (0..$nchunk)
{

  print "chunk# $chunk \n";

  print "Press enter to continue \n";
  my $inp = <>;

  for my $work (0..$chunksize -1)
  {
    print "work# $work of $chunksize \n";

    my $idx = $chunk * $chunksize + $work ; 
    if( $idx < $#kys)
    {
      my $stem = $kys[$idx]; 
      my $f =   "tmp/${stem}_plot.gp";

      if( -e $f )
      {
        unlink $f; 
      }

      open GNU , ">" , $f;  
      print GNU "set multiplot layout 2 , 2 rowsfirst title \"";

      my @chars = map substr($stem, $_, 1), 0 .. length($stem) -1; 
      my $idx = 0;
      foreach my $char (@chars) 
      { 
        if($idx == 80)
        {
          print GNU "\\n ";
          $idx = 0; 
        }
        print GNU $char;
        $idx++; 
      }
      print GNU "\" \n";


      my $r1 = gen_dat($stem,".dat_corr_pre");  
      my @pre = @{ $r1 };
      print_gnu("pre" , @pre) ; 

      my $r2 = gen_dat($stem,".dat_corr_post");  
      my @post = @{ $r2 };
      print_gnu("post" , @post); 

      my $r3 = gen_dat($stem,".dat_norm");  
      my @norm = @{ $r3 };
      print_gnu("norm",@norm);  

      print GNU "unset multiplot \n";
      close GNU; 

      system ("gnuplot -persist $f "); 
    }
  }

}



sub gen_dat
{
  my ($stem,$suffix) = @_; 

  my $f = $stem . $suffix; 

  my $r = "tmp/${f}.real" ;
  my $i =  "tmp/${f}.imag"; 


  if( -e $r )
  {
    unlink $r; 
  }

  if( -e $i )
  {
    unlink $i;
  }

  my $cmdr = " calcbc \' real ( $f ) \' > $r "; 
  my $cmdi = " calcbc \' imag ( $f ) \' > $i ";

  system( $cmdr ) == 0 || die($_); 
  system( $cmdi ) == 0 || die($_); 

  my @fs = (); 
  push @fs , $r; 
  push @fs , $i;

  return \@fs; 
}


sub print_gnu
{
  my ($ti,$rl,$im) = @_; 
  print GNU "set label 1 \"$ti\" at graph 0.45,0.9 font \',8\' \n";
  print GNU "plot \'$rl\' using 1:2:3 w yerr ti \'real\' , \\\n";
  print GNU "\'$im\' using 1:2:3 w yerr ti \' imag\' \n";
}
