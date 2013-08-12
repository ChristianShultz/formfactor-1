#!/bin/tcsh

source $HOME/.MPI_ENV.csh
source $HOME/.CUDA_ENV.csh

set xml="dist_src.ini.xml"
set src="src.mod"

/bin/rm -f $src ${xml}.out


if( ! -f $xml ) then
  echo " distillution xml ini is not present : $xml "
  exit 1
endif

set harom=$HOME/git-builds/harom-jitptx/gpu-parscalar-Nd3/bin/harom

if( ! -e $harom ) then 
  echo " harom does not exist : $harom " 
  exit 1
endif 

$harom -i $xml >&! ${xml}.out

if ( $status != 0 ) then 
  echo " harom distillution failed "
  exit 1
endif 

