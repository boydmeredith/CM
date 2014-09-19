#!/bin/bash
#$ -M jtboyd@stanford.edu
#$ -m beas
#$ -N cm_job

for trWs in "[0 0 .33 .34 .33 0]" "[1 0 0 0 0 0 ]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]"  "[0 0 0 1 0 0]"  "[0 0 0 0 1 0]" "[0 0 0 0 0 1]" 
do
 

for which_trte in  "[1 2 3 4 0 0 0 0]" "[0 0 0 0 1 2 3 4]" "[1 1 1 2 2 2 2 2]"
do

for cNames in "{'hit', 'cr'}" "{'s_old','s_new'}" "{'hit', 'fa'}" "{'cr', 'miss'}" "{'old',new'}" "{'fa','miss'}" "{'cr', 'fa'}" "{'hit', 'miss'}" "{'hit','miss','fa','cr'}"
do

qsub -l h_vmem=10G -v -cwd cm_job.bash "[],'','trWeights',${trWs},'condNames',${cNames},'which_traintest',${which_trte}" 


done
done
done
