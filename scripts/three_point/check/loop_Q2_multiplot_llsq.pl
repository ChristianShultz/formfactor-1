#!/usr/bin/perl -w


BEGIN{
  push (@INC, "/u/home/shultz/git/radmat/scripts/three_point/check");
}

use strict; 
use findPhase; 

# assumes being run from the top

my @qs = split (/\s+/,  ` ls -al | grep ^d | awk '{ print \$9 }' | grep Q2_ | xargs ` );


my $base = `pwd`;

chomp $base; 

foreach my $q2 (@qs)
{
  chomp $q2;

  print "type return to continue \n"; 

  my $foo = <>; 

  chdir $base || die (" couldnt move to $base " ); 

  # do things that are not zero

  my $q2dir =  $base . "/" . $q2;
  chdir $q2dir || die(" coudnt move to $q2 "); 


  my @lgs = split (/\s+/,  ` ls -al | grep ^d | awk '{ print \$9 }' | grep lefty | xargs ` );

  foreach my $lg (@lgs)
  {
    chdir $q2dir || die(" coudnt move to $q2 "); 

    my $lgdir =   $lg . "/llsq" ;
    print $lgdir;

    if ( ! -d $lgdir )
    {
      print " warning $lg has no llsq dir. \n";
      next; 
    }
    chdir $lgdir || die(" coudnt move to $lg "); 

    my $tag = `head -1 row_index_to_continuum_elem.txt | awk \' { print \$2 } \' `;
    print $tag ."\n";
    my $ph = findPhase->new(); 
    $ph->id($tag); 
    $ph->parse_id();
    print "pf = " . $ph->psnk_str() . "\nq = " . $ph->q_str() . "\npi = " . $ph->psrc_str() . "\n"; 



    print "Q2 = $q2 things that should be non zero \n"; 

    my @tst = split (/\s+/, `ls *.jack | xargs` ); 

    if( @tst )
    {
      system ("multiplot_all_jack_inp.pl 1 1") == 0 || die (" multiplot failed "); 
    }
    else
    {
      print "directory was empty, skipping \n";
    }

    print "type return to continue \n";

    my $foofoo = <>;

    my $zero_dir = "zeroed_matrix_elems";

    chdir $zero_dir || die(" couldnt move to $zero_dir \n" ); 

    print "Q2 = $q2 things that should be zero \n"; 

    @tst = split (/\s+/, `ls *.jack | xargs` ); 

    if( @tst )
    {
      system ("multiplot_all_jack_inp.pl 1 1") == 0 || die (" multiplot failed "); 
    }
    else
    {
      print "directory was empty, skipping \n";
    }
  }
}



