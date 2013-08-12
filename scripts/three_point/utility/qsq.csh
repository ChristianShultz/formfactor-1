#!/bin/csh

if($#argv != 10) then 
echo "$0 : usage : m_f m_i pfx pfy pfz pix piy piz ls xi "
exit 1
endif

set mf=$1
set mi=$2
set pfx=$3
set pfy=$4
set pfz=$5
set pix=$6
set piy=$7
set piz=$8
set ls=$9
set xi=$10
set tp="6.2831"

# get pfsq
set cmd1=" $pfx * $pfx + $pfy * $pfy + $pfz * $pfz ";
set ppf=`eval echo " '$cmd1' | bc -l "`

# get pisq in
set cmd2=" $pix * $pix + $piy * $piy + $piz * $piz ";
set ppi=`eval echo " '$cmd2' | bc -l "`

#get sapce q
set cmd3=" ( ($pfx - $pix)^2  + ($pfy - $piy)^2  + ($pfz - $piz)^2 ) *  ( $tp / $xi / $ls )^2    " 
set q2space=`eval echo " '$cmd3' | bc -l "`

set Ef=`boost.csh $mf $ppf $ls $xi`
set Ei=`boost.csh $mi $ppi $ls $xi`

# get temporal q
set cmd4=" ($Ef - $Ei)^2 "
set q2time=`eval echo " '$cmd4' | bc -l "`

set cmd5=" $q2space - $q2time ";
set q2=`eval echo " '$cmd5' | bc -l "`

echo $q2

