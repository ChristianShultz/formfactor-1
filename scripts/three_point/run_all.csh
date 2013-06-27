#!/bin/tcsh


source $HOME/.MPI_ENV.csh
source $HOME/.CUDA_ENV.csh

# clean up the rubbish 
# ./clean.csh 

# make a distillution source
# ./run_dist_src_harom-jitptx.csh 

# generate solution vectors 
# ./run_prop_dist_chroma-jitptx.csh 




# actually do something 

foreach x ( clean.csh \
            run_dist_src_harom-jitptx.csh \
            run_prop_dist_chroma-jitptx.csh \
            run_peram_dist_harom-jitptx.csh \
          )

       ./${x}

      if ($status != 0) then 
        echo " $x failed " 
        exit 1
      endif

end


