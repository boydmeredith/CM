#!/bin/bash

#$ -N cm_job

for trWeights in "[1 0 0 0 0 0]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]"  "[0 0 0 1 0 0]"  "[0 0 0 0 1 0]" "[0 0 0 0 0 1]" "[0 0 .33 .34 .33 0]"
do
for condNames in "{'s_old',s_new'}" "{'old','new'}" "{'hit', 'cr'}" "{'fa','miss'}" "{'cr', 'fa'}" "{'cr', 'miss'}" "{'hit', 'miss'}" "{'hit', 'fa'}" "{'hit','miss','fa','cr'}"
do
for which_traintest in "[1 2 3 4 0 0 0 0]" "[0 0 0 0 1 2 3 4]" "[1 2 3 4 1 2 3 4]" "[1 1 2 2 3 3 4 4]" "[1 2 1 2 3 4 3 4]" "[1 1 1 1 2 2 2 2]"
do
qsub -l h_vmem=10G -v trWeights=$trWeights condNames=$condNames which_traintest=$which_traintest -cwd cm_job.bash
done
done
done
