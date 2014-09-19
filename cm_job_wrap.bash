#!/bin/bash
#$ -M jtboyd@stanford.edu
#$ -m beas
#$ -N cm_job

for cNames in "{'old','new'}" "{'hit', 'fa'}" "{'s_old',s_new'}" "{'hit', 'cr'}" "{'fa','miss'}" "{'cr', 'fa'}" "{'cr', 'miss'}" "{'hit', 'miss'}" "{'hit','miss','fa','cr'}"
do

for trWs in "[0 0 .33 .34 .33 0]" #"[1 0 0 0 0 0]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]"  "[0 0 0 1 0 0]"  "[0 0 0 0 1 0]" "[0 0 0 0 0 1]" 
do

for which_trte in "[1 2 3 4 5 6 7 8]" #"[1 2 3 4 0 0 0 0]" "[0 0 0 0 1 2 3 4]" "[1 2 3 4 1 2 3 4]" "[1 1 2 2 3 3 4 4]" "[1 2 1 2 3 4 3 4]" "[1 1 1 1 2 2 2 2]"
do

 

qsub -l h_vmem=10G -v -cwd cm_job.bash "[],'','trWeights',${trWs},'condNames',${cNames},'which_traintest',${which_trte}" 


done
done
done
