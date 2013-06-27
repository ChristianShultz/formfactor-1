#!/bin/tcsh

source $HOME/.MPI_ENV.csh
source $HOME/.CUDA_ENV.csh


set xml="prop_dist.ini.xml"
set prop="prop.mod"

/bin/rm  -f $prop ${xml}.out

if (! -f $xml ) then 
  echo " prop distillution ini is not present : $xml "
  exit 1
endif

set chroma=$HOME/git-builds/chroma-jitptx/parscalar-jit-Nd4/bin/chroma

if( ! -e $chroma ) then 
  echo " chroma does not exist : $chroma "
  exit 1
endif

$chroma -i $xml >&! ${xml}.out

if ( $status != 0 ) then 
  echo " chroma prop dist failed "
  exit 1
endif


