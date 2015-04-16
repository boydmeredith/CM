
module load MATLAB-R2012b
module avail
echo "CM_run_mvpa_srch($1,'$2',1, 'trWeights_train', $3, 'trWeights_test', $4)"
matlab -nodesktop -r "CM_run_mvpa_srch($1,'$2',1, 'trWeights_train', $3, 'trWeights_test', $4)"
