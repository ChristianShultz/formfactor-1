#!/bin/tcsh

source $HOME/.MPI_ENV.csh
source $HOME/.CUDA_ENV.csh


/bin/rm -f hadron_node.sdb1

/home/shultz/git-builds/harom-jitptx/gpu-parscalar-Nd3/bin/harom -i harom_hadron_node.xml
