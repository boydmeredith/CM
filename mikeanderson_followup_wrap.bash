qsub -l h_vmem=10G -v -cwd mikeanderson_followup.bash '[1 2 2 2 1 1 1 1]' 1
qsub -l h_vmem=10G -v -cwd mikeanderson_followup.bash '[2 2 2 1 1 1 1 1]' 1

