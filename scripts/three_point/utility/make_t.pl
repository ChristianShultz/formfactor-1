#!/usr/bin/perl

die "make_t.pl <cfgs> <Lt> \n" unless $#ARGV == 1; 
$cfgs = $ARGV[0];
$Lt = $ARGV[1];


open(T, "> t");
print T "$cfgs $Lt 0 0 1\n";
for($i=0; $i < $cfgs; $i++){
    for($t=0; $t < $Lt; $t++){
        print T "$t $t\n";
    }
}
close(T);

exit(0);
