module load MATLAB-2012b
module avail
matlab -nodesktop -r "CM_run_mvpa_v4_tbm([], '', 'which_traintest', $1, 'condNames', {'exHit', 'exCr'}, 'scramble', $2)"

