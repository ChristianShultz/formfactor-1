#!/usr/bin/perl -w

BEGIN{
  push (@INC, "/Users/shultz/jlab/scripts");
}



use strict;
use Math::Complex;
use Math::Complex ':pi';
use findPhase; 
use threads;
use threads::shared;




die("usage: $0 <MatrixElement_> <KinematicZeroCut>") unless $#ARGV == 1; 

my $matrix_elememt_stem = $ARGV[0];
my $kcut = $ARGV[1]; 

#
# what is numerically zero ?
my $tol_deg = 15;
my $tol = $tol_deg * 3.1415 / 180.; 
my $tol_per_pi = $tol / 3.1415;


print " using tolerance $tol = $tol_deg (deg) = $tol_per_pi * pi\n"; 

#
# where do we want to look for the phase
my $tlow = 8;
my $thigh = 12; 


# assumption that you are running in the llsq directory
my $llsq_file = "row_index_to_continuum_elem.txt";

my $foo = `ls -al | grep ^d  | grep Q2 | grep -v FF | awk \'{print \$9} \' | xargs `;

my @QS = split /\s+/ , $foo ; 

foreach my $q (@QS) 
{
  $q = $q . "/llsq";
}

my @bad_elems;  

my $dir = `pwd`;
chomp $dir;

foreach my $i (0..$#QS)
{
  chdir $dir || die("couldn't move back up");
  my $ref = find_phase($QS[$i]);
  my @a = @{ $ref };

  if(@a)
  {
    print "Q2 = $QS[$i]  had $#a baddies \n";

    push @bad_elems , @a; 

    print " bad phases \n" ;

    foreach my $elem (@a)
    {
      print $elem->elem_out_str() . "\n" ;
    }

  }  
}


if (@bad_elems)
{
  print " had $#bad_elems bad phases \n" ;
}
else
{
  print " all phases passed \n";
}


###########################################################


sub find_phase
{

  my $dir = shift; 

  chdir $dir; 

  print "working on $dir \n"; 

# init the array of findPhase "classes"
  my @elems = (); 

  return \@elems unless -e $llsq_file;

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

  my @kcut_elems = @elems; 
  @elems = (); 
  @elems = @{ run_kinematic_cut($kcut,@kcut_elems) };


# loop and check the phases
  foreach my $elem (@elems)
  {

    my $zlat = Math::Complex->make($elem->lat_real(),$elem->lat_imag()); 
    my $zrad = Math::Complex->make($elem->rad_real(),$elem->rad_imag());
    my $philat;
    my $phirad;

      $philat = sprintf( "%.3f",  arg($zlat) - $philat0 ); 
      $phirad = sprintf( "%.3f" , arg($zrad) - $phirad0 ); 

    $elem->philat($philat);
    $elem->phirad($phirad); 

  }


  my @baddies = (); 

  foreach my $elem (@elems)
  {
    if( abs ( $elem->phase_diff() ) > $tol )
    {
      # kill 2pi
      if( abs(abs ( $elem->phase_diff() ) - 2. * 3.1415) > $tol)
      {
        push @baddies, $elem; 
      }
    }
  }

  return \@baddies; 
}


sub run_kinematic_cut
{

  my ($k,@eles) = @_; 
  my @elems = (); 

  foreach my $elem (@eles)
  {

    my $zlat = Math::Complex->make($elem->lat_real(),$elem->lat_imag()); 
    my $zrad = Math::Complex->make($elem->rad_real(),$elem->rad_imag());
    
    my $rlat = abs($zlat); 
    my $rrad = abs($zrad); 

    if( ($rlat > $k) && ($rrad > $k) ) 
    {
      push @elems , $elem; 
    }
    elsif ($rlat > $k) 
    {
      print " FOUND A STRANGE KCUT norm(lat) = $rlat norm(rad) = $rrad \n";
    } 
    elsif ($rrad > $k) 
    {
      print " FOUND A STRANGE KCUT norm(lat) = $rlat norm(rad) = $rrad \n";
    }
    else
    {
      print " KCUT norm(lat) = $rlat norm(rad) = $rrad cut = $k \n";
    } 
  }

  return \@elems; 
}



