#$ -cwd
#!/bin/bash

#$ -M $jtboyd@stanford.edu
#$ -m beas

#$ -N cm_mvpa

module load MATLAB-R2012b
echo ${foo}

matlab -r "CM_run_mvpa_v4_tbm(${foo})"

