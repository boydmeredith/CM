for roiName in ""
do
for trte in "[1 2 3 4 0 0 0 0]" "[1 1 1 1 2 2 2 2]" "[0 0 0 0 1 2 3 4]" "[4 3 2 1 1 2 3 4]" "[1 2 3 4 1 2 3 4]"
do	
for trW in "[0 0 0 0 0 1]" "[0 0 0 0 1 0]" "[0 0 0 1 0 0]" "[0 0 1 0 0 0]" "[0 1 0 0 0 0]" "[1 0 0 0 0 0]"
do
/Applications/MATLAB_R2011b.app/bin/matlab -nodesktop -r "CM_run_mvpa_v4_tbm([],'','roiName',${roiName},'trWeights',${trW}, 'which_traintest', ${trte})"
/Applications/MATLAB_R2011b.app/bin/matlab -nodesktop -r "CM_run_mvpa_v4_tbm([],'','roiName',${roiName},'trWeights',${trW}, 'which_traintest', ${trte},'condNames',{'leftHand','rightHand'})"
done
