
module load MATLAB-2012b
module avail
echo "CM_run_mvpa_srch($1,'$3',1)"
matlab -nodesktop -r "CM_run_mvpa_srch($1,'$3',1)"
