#!/bin/tcsh -x

if ($#argv != 2 ) then 
  echo "usage: <file> <nsub>"
  exit 1
endif

set f=$argv[1]
set y=$argv[2]

foreach x ( `seq 1 1 $y` )
  qsub $f 
  sleep 1
end


