#!/bin/tcsh -x


if ($#argv != 1) then 
  echo "usage: $0 <dt> "
  exit 1
endif

set dt = $argv[1]

# run radmat_util

set rad="unsmeared._DT_.ini.xml"
set radt="unsmeared.${dt}.ini.xml"

if ( ! -f $rad ) then 
  echo "missing $rad"
  exit 1
endif


# replace with dt
cat $rad | awk -v TIME=$dt '{ gsub( "_DT_" , TIME ); print}' > $radt
set rade=/home/shultz/git-builds/radmat/bin/radmat_util

if (! -e $rade ) then 
  echo "missing $rade"
  exit 1
endif

# generate the npt lists
$rade gen_xml $radt split_mom 

# strip the first line which is some <?xml nonsense>
foreach xml ( `ls | grep "npt.list"` )
  
  tail -n +2 $xml > tmp.xml 
  mv tmp.xml $xml 

end

set cfgs="sdb.list"

if (! -f $cfgs) then 
  echo "missing $cfgs"
  exit 1
endif

set skel="ib.8c8.skel.csh"

foreach x ( p000 \
            p100 \
            p110 \
            p111 \
            p200 \
          )

  echo $x

  # create a run script 
  set run="ib.8c8.${x}.csh"
  cat $skel | awk -v PP=$x -v TIME=$dt '{ gsub ("_DT_",TIME) ; gsub("_MOM_",PP); print}' > $run
  chmod 0755 $run


  # make a run dir for this momentum 
  mkdir $x

  # copy the sdblist
  cp $cfgs $x/$cfgs
  
  # my momentum npt
  set npt="npt.list.${x}.xml"; 

  if (! -f $npt) then 
    echo "missing $npt"
    exit 1
  endif

  # make a local copy of it
  cp $npt $x/$npt

  # move down intp  pXYZ
  pushd $x

  # link the relevant programs so we can just do a local run 
  ln -s ../ParamsDistillation.pm
  ln -s ../drive_unsmeared_nodes.pXYZ.pl
  ln -s ../submit.csh
  ln -s ../$run

  # submit 60 jobs for this momentum (they auto resub until the sequence is finished)
  ./submit.csh $run 60 & 

  # jump back up
  popd

end

