thresh=.01
for i in {2..6}; do
analysis=trByTr_exToCm_hitVCrs_wholebrain_crossTr_tr${i}ToTr${i}
echo $analysis 
qsub -l h_vmem=10G -v -cwd srch_vis_job.bash $analysis $thresh
done
   
