for i in 1 `seq 3 10` `seq 12 26`; do
qsub -l h_vmem=10G -v -cwd cm_sl_mvpa_single_sub.bash $i $1
done
 
