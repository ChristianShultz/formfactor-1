#!/bin/tcsh -x

setenv PBS_NODEFILE "nodefile.node"

if ( -f $PBS_NODEFILE ) then 
  /bin/rm $PBS_NODEFILE 
endif

echo "qcd12kmi" > $PBS_NODEFILE

source $HOME/.MPI_ENV.csh

set seqno=1000 
set dt=28 
set outp=$HOME/shultz_volatile/szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per
#set outp=.

./clean_scratch.sh

 ./drive_unsmeared_nodes.pXYZ.pl $seqno $dt $outp p000 
