#!/usr/bin/perl -w

BEGIN{
  push (@INC, "/Users/shultz/jlab/scripts");
}



use strict;
use Math::Complex;
use Math::Complex ':pi';
use findPhase; 


#
# sanity 
#
# my $z = Math::Complex->make(5,6);
# should be ~0.876 in rads ~50 in deg
# print arg($z) . " \n" ;


# my $z = Math::Complex->make(0.,0.);
# print arg($z) . "\n";
# die("unholy death");

die("usage: $0 <MatrixElement_> <tlow> <thigh>") unless $#ARGV == 2; 

my $matrix_elememt_stem = $ARGV[0];

#
# what is numerically zero ?
my $toldeg = 15.; 
my $tol = $toldeg*3.1415/180.; 


print "using tolerance $tol ~ $toldeg (deg) \n";
#
# where do we want to look for the phase
my $tlow = $ARGV[1];
my $thigh = $ARGV[2]; 


# assumption that you are running in the llsq directory
my $llsq_file = "row_index_to_continuum_elem.txt";

# init the array of findPhase "classes"
my @elems = (); 

open IN , "<" , $llsq_file;
while (<IN>)
{
  my ($elem , $id) = split(/\s+/, $_);  
  my $phase = findPhase->new(); 
  $phase->matrix_element_stem( $matrix_elememt_stem );
  $phase->elem($elem); 
  $phase->id($id);
  $phase->parse_id();  
  push @elems , $phase; 
}


# fill in the bit calculated by radmat
foreach my $elem (@elems)
{
  my $cmd = "grab_ff " . $elem->rad_str();
  my $opt = `  $cmd  `;
  my ($rl, $im) = split(/\s+/,$opt); 
  $elem->rad_real($rl); 
  $elem->rad_imag($im); 
}

#  die("unholy death");

# fill in the real and imag parts from the lattice 
foreach my $elem (@elems)
{
  my $cmd = "grab_vector_mean_cmp " . $elem->lat_elem();
  # print $cmd . "\n";
  my $opt = `  $cmd $tlow $thigh `;
  # print $opt . "\n";

  my ($rl , $im ) = split(/\s+/,$opt);
  $elem->lat_real($rl);
  $elem->lat_imag($im); 
}


# shift an element to get a reference 
my $elem0 = shift( @elems ) ; 
my $zlat0 = Math::Complex->make($elem0->lat_real(),$elem0->lat_imag()); 
my $zrad0 = Math::Complex->make($elem0->rad_real(),$elem0->rad_imag());
my $philat0 = arg($zlat0);
my $phirad0 = arg($zrad0); 

# push it back into the list since I need closure
push @elems , $elem0; 

# loop and check the phases
foreach my $elem (@elems)
{

  my $zlat = Math::Complex->make($elem->lat_real(),$elem->lat_imag()); 
  my $zrad = Math::Complex->make($elem->rad_real(),$elem->rad_imag());

#  phase difference  
  my $philat = sprintf( "%.3f",  arg($zlat) - $philat0 ); 
  my $phirad = sprintf( "%.3f" , arg($zrad) - $phirad0 ); 

#  full phase
#  my $philat = sprintf( "%.3f",  arg($zlat) ); 
#  my $phirad = sprintf( "%.3f" , arg($zrad) ); 


  $elem->philat($philat);
  $elem->phirad($phirad); 
}


# screen dump

@elems = sort {$a->elem() <=> $b->elem} @elems; 

# @elems = sort {$a->psnk() < $b->psnk()} @elems;

my @superbad = (); 

foreach my $elem (@elems)
{

  if("true")
  {
    if( abs ( $elem->phase_diff() ) > $tol )
    {
      # kill 2pi
      if( abs ( abs ( $elem->phase_diff() ) - 2. * 3.1415)  > $tol )
      {
        #    by elem 
        print "*" . $elem->elem_out_str() . "\n";  

        push @superbad , $elem; 
        #    by tag 
        #   print $elem->out_str() . "\n";
      }
      else
      {
        print " " . $elem->elem_out_str() . "\n";  
      }

    }
    else
    {
      print " " . $elem->elem_out_str() . "\n";  
    }
  }
  else
  {
    print $elem->elem_out_str() . "\n";  
  }
}


print "******** \nbad list \n********\n";
foreach my $elem (@superbad)
{
        print "*" . $elem->elem_out_str() . "\n";  
}
