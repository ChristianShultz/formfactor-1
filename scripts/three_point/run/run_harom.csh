#!/bin/tcsh -x



if ($#argv != 1) then
  echo "Usage: $0  <redstar xml input file>"
  exit 1
endif

set t="harom_hadron_node.xml"
set f=$1

if (! -f $t) then
  echo "meson_hadron_node template file does not exist: $t"
  exit 1
endif

if (! -f $f) then
  echo "redstar input file does not exist: $f"
  exit 1
endif

set g=${f:r}

/bin/rm -f hadron_node.xml

source ~/.MPI_ENV.csh
source ~/.CUDA_ENV.csh


set exe_red=/home/shultz/git-builds/redstar/bin
set exe_harom=/home/shultz/git-builds/harom/scalar-Nd3/bin


$exe_red/redstar_gen_graph $f out.xml >&! ${g}.gen.out
if ($status != 0) then
  echo "redstar_gen_graph failed"
  exit 1
endif

$exe_harom/harom -i $t -o out.xml >&! had.${f}.out
if ($status != 0) then
  echo "harom failed"
  exit 1
endif

# removed from orig.. don't know what tmp was suppose to be
##/bin/rm $tmp

$exe_red/redstar_npt $f out.xml >&! ${g}.out
if ($status != 0) then
  echo "redstar_npt failed"
  exit 1
endif


