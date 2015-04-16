module load MATLAB-R2012b
module avail

matlab -nodesktop -r "sl_mvpa_visualize('$1', $2, 1)"
