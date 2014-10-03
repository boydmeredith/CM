#!/bin/bash

#$ -N matlab_CMMVPA

module load MATLAB-R2012b
matlab -r "cd ~/cm/Data/Functional/mvpa_results; CM_batchMvpaPostProc(pwd, 'conds*1_2_3_4*',[],'')"
