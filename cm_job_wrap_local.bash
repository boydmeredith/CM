#!/bin/bash
#$ -M jtboyd@stanford.edu
#$ -m beas
#$ -N cm_job

for cNames in "{'hit','cr'}" #"{'hit','fa'}" "{'cr','miss'}" #"{'fa','miss'}" "{'cr', 'fa'}" "{'cr', 'miss'}" "{'hit', 'miss'}" "{'hit','miss','fa','cr'}" "{'old','new'}" 
do

for trWs in "[0 0 .33 .34 .33 0]" "[1 0 0 0 0 0]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]"  "[0 0 0 1 0 0]"  "[0 0 0 0 1 0]" "[0 0 0 0 0 1]" 
do

for which_trte in "[1 1 1 1 2 2 2 2]" #"[1 2 3 4 1 2 3 4]" "[1 1 2 2 3 3 4 4]" "[1 2 1 2 3 4 3 4]" "[1 1 1 1 2 2 2 2]"
do

 

	
/Applications/MATLAB_R2011b.app/bin/matlab -nodesktop -r  "CM_run_mvpa_v4_tbm([],'','trWeights',${trWs},'condNames',${cNames},'which_traintest',${which_trte})" 
/Applications/MATLAB_R2011b.app/bin/matlab -nodesktop -r "exit"

done
done
done
