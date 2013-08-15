#!/bin/tcsh -x
#PBS -N dt28u_p200_c8
#PBS -q ib
#PBS -A Spectrum
#PBS -l walltime=48:00:00
#PBS -l nodes=8:cores8
#PBS -j eo


source $HOME/.JITPTX_ENV.csh

if ($?PBS_O_HOST != 0) then
  cd $PBS_O_WORKDIR
  echo "PBS_JOBID = $PBS_JOBID"
  echo "PBS_JOBNAME = $PBS_JOBNAME" 
  echo "PBS_O_WORKDIR = $PBS_O_WORKDIR"
  echo "PBS_NODEFILE"
  cat $PBS_NODEFILE | uniq
endif

set file_prefix="szscl3_16_128_b1p50_t_x4p300_um0p0743_n1p265_per"
set list="sdb.list"
set run_script="${PBS_O_WORKDIR}/drive_unsmeared_nodes.pXYZ.pl"
set dt=28
set mom="p200"

if ( ! -e $run_script ) then
  echo "$run_script is not present"
  exit 0
endif


if (! -d "/volatile/Spectrum") then
  echo "Lustre is not present"
  exit 0
endif

set chunk=1

set nn=0

while ($nn < $chunk)

  #
  # Decide if the job needs to be resubmitted.
  #
  if ( ! -f ${list} ) then
    echo "ERROR: file  $list  does not exist"
    exit(11)
  endif

  set numseqno=`cat $list | wc -l`

  if ( $numseqno == 0 ) then
    touch ${file_prefix}.done
    echo "The $file_prefix sequence is finished"
    exit(0)
  endif

  set newseqno=`head -1 $list | gawk '{print $NF}'`
  #set newseqno = 1040

  # Bump up to the next item in the list
  cat $list | awk 'NR>1{print $0}' > list_$$
  /bin/mv -f list_$$ $list


  set output_path = $HOME/shultz_volatile/$file_prefix;
  set out = $output_path/${file_prefix}.out.dt${dt}.${mom}.${newseqno}
  /bin/rm -f $out ${out}.gz
  touch $out

  echo "Stem= $file_prefix  Cfg = $newseqno"
  echo "Stem= $file_prefix  Cfg = $newseqno" >> $out

  if ($?PBS_O_HOST != 0) then
    echo "PBS_JOBID = $PBS_JOBID" >> $out
    echo "PBS_JOBNAME = $PBS_JOBNAME"  >> $out
    echo "PBS_O_WORKDIR = $PBS_O_WORKDIR" >> $out
    echo "PBS_NODEFILE" >> $out
    cat $PBS_NODEFILE | uniq >> $out
  endif

  echo "DELTAT = $dt "  
  echo "DELTAT = $dt " >> $out 

  echo "/scratch/ before run " >> $out
  ls -al /scratch >> $out

  time $run_script $newseqno $dt $output_path $mom >> $out

  if ( $status != 0) then
    echo "Some error in run for $out at seqno=$newseqno"
    echo $newseqno >> ${list}.badlist
    exit 1
  endif
  echo "Ending chroma script at time= " `date`
  echo "Ending chroma script at time= " `date` >> $out


  echo "/scratch/ post run " >> $out
  ls -al /scratch >> $out

  # using jit these grow to ~25mb
  # /bin/rm -f $out ${out}.gz
  # gzip -f9 $out

  @ nn++
end   # while nn

if ($?PBS_O_WORKDIR != 0) then
  /usr/local/bin/qsub ib.8c8.p200.csh
endif

exit 0

