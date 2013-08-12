#!/bin/bash

export LD_LIBRARY_PATH=/dist/gcc-4.6.3/lib64:/dist/gcc-4.6.3/lib:$LD_LIBRARY_PATH

ARGS=$*

case $MV2_COMM_WORLD_LOCAL_RANK in
0)
 CPUS="0,1"
 ;;
1)
 CPUS="2,3"
 ;; 
2)
 CPUS="4,5"
 ;;
3)
 CPUS="6,7"
 ;;
*) 
  echo LOCAL Rank cannot be bigger than 3
  exit 1
  ;;
esac

/usr/bin/numactl --physcpubind=${CPUS} ${ARGS}
