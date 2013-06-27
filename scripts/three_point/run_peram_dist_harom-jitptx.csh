#!/bin/tcsh

source $HOME/.MPI_ENV.csh
source $HOME/.CUDA_ENV.csh

set xml="peram_dist.ini.xml"
set src="peram.mod"

/bin/rm -f $src ${xml}.out


if( ! -f $xml ) then
  echo " peram distillution xml ini is not present : $xml "
  exit 1
endif

set harom=$HOME/git-builds/harom-jitptx/gpu-parscalar-Nd3/bin/harom

if( ! -e $harom ) then 
  echo " harom does not exist : $harom " 
  exit 1
endif 

$harom -i $xml >&! ${xml}.out

if ( $status != 0 ) then 
  echo " harom peram distillution failed "
  exit 1
endif 

