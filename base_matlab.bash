#!/bin/bash

#$ -N matlab_CMMVPA

module load MATLAB-R2012b
matlab -r "CM_run_mvpa_v4_tbm([],'leftIndexRightMiddle_00112222_motorMask','trWeights',${foo},'condNames', {'leftIndex','rightMiddle'}, 'which_traintest', [0 0 1 1 2 2 2 2], 'roiName', 'rLR_Precentral_aal')"
