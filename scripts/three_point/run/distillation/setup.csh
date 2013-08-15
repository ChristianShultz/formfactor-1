#!/bin/tcsh -x

set cfgs="sdb.list"

if (! -f $cfgs) then 
  echo "missing $cfgs"
  exit 1
endif

foreach x ( p000 \
            p100 \
            p110 \
            p111 \
            p200 \
          )

  echo $x

  mkdir $x
  cp $cfgs $x/$cfgs
  
  set npt="npt.list.${x}.xml"; 

  if (! -f $npt) then 
    echo "missing $npt"
    exit 1
  endif

  cp $npt $x/$npt

  set run="ib.8c8.${x}.csh"

  if (! -f $run) then 
    echo "missing $run"
    exit 1
  endif

  pushd $x

  ln -s ../ParamsDistillation.pm
  ln -s ../drive_unsmeared_nodes.pXYZ.pl
  ln -s ../submit.csh
  ln -s ../$run

  ./submit.csh $run 60 & 

  popd

end

