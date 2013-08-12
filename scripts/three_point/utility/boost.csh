#!/bin/csh 

if($#argv != 4) then 
  echo "$0 : usage : m_rest pp ls xi "
  exit 1
endif


set m=$1;
set pp=$2;
set ls=$3
set xi=$4
set tp="6.2831"

set cmd=" sqrt( $m * $m +  ( $tp / $xi / $ls )^2 * $pp ) "
eval echo " '$cmd'  | bc -l"
