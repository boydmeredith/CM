for tr_train in "[1 0 0 0 0 0]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]" "[0 0 0 1 0 0]" "[0 0 0 0 1 0]" "[0 0 0 0 0 1]"; do
for tr_test in "[1 0 0 0 0 0]" "[0 1 0 0 0 0]" "[0 0 1 0 0 0]" "[0 0 0 1 0 0]" "[0 0 0 0 1 0]" "[0 0 0 0 0 1]"; do
for i in 1 `seq 3 10` `seq 12 26`; do
qsub -l h_vmem=10G -v -cwd cm_sl_mvpa_single_sub.bash $i $1 "$tr_train" "$tr_test"
done
done
done
