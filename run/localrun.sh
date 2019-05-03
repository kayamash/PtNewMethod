#!/bin/sh
run="nohup root -l -b -q run.cpp++ > /home/kayamash/log/LUT.txt &"
clean="rm -r run_*"
eval $run
eval $clean
