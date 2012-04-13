#!/usr/bin/perl

# test_all.pl -
#
# Tuesday, April  3 2012
#


 @exes = qx (find * -type f -perm -o+rx);

pop @exes;


foreach (@exes)
  {
    print "testing $_ ";
    system ("./$_") == 0 || die "$_ failed";
    print "\n\n";
  }

#system (@exes) == 0 || die "system @exes failed: $?";
