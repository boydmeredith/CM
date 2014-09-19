#!/bin/bash

for i in "[1 0 0 0 0 0]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]" "[0 0 0 1 0 0]" "[0 0 0 0 1 0]" "[0 0 0 0 0 1]"
do

qsub -l h_vmem=10G -v foo=$i /hsgs/projects/awagner/jtboyd/CM_ret/Scripts  base_matlab.bash

done
