#!/usr/bin/perl
use strict;

die "forall.pl <t0> <tZ> <state> <opslistfile> <opname> \n" unless $#ARGV == 4;

my $t0 = $ARGV[0];
my $tZ = $ARGV[1];
my $state = $ARGV[2];
my $opslistfile = $ARGV[3];
my $opname = $ARGV[4];

my $massfile = "t0${t0}/MassJackFiles/mass_t0_${t0}_reorder_state${state}.jack";

my $ii = 0;
open(OPSLIST, "< $opslistfile") || die("Could not open ${opslistfile}\n");
while(my $subopnameline = <OPSLIST>)
{
  chomp $subopnameline;
  my @subopnameline_split = split(/ /, "${subopnameline}");
  my $subopname = $subopnameline_split[1];

  my $file = "t0${t0}/V_tJackFiles/V_t0_${t0}_reordered_state${state}_op${ii}.jack";

  my $command = "calcbc \" sqrt( 2 * ${massfile} ) * exp ( - ${massfile} * ${t0} / 2 ) * extract ( ${file} , ${tZ} ) \" | awk '{print \$2}' ";
  my $val = `${command}`;
  print "${opname} ${subopname} ${val}";
  $ii++;
}
close(OPSLIST);

