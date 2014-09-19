module load MATLAB-R2012b
for trW in "[0 1 0 0 0 0]"
do
foo="[],'','trWeights',${trW}"
#bash cm_job.bash ${foo}
qsub -l h_vmem=10G -v foo="[],'','trWeihts',${trW}" -cwd cm_job.bash "CM_run_mvpa_v4_tbm([],'','trWeights',${trW})"
done
