#!/bin/bash
run="root -l -b -q run.cpp++"
bsub="bsub -q 4h -o ~/PtNewMethodLog/out.log -e ~/PtNewMethodLog/err.log "
clean="rm run_*"
reset="rm *.log"
eval $reset
eval $bsub$run
eval $clean
