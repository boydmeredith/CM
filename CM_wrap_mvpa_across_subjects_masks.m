function [res] = CM_wrap_mvpa_across_subjects_masks()
%CM_wrap_mvpa_across_subjects_masks() is intended to run classification
%across subjects using pre-specified masks

subj_array = [1 3:10 12:26];

for s = subj_array
    	
	trte_str = sprintf('minus%02d',s);
	mask_to_use = sprintf('exhitsgrcrs_AND_cmhitsgrcrs_minus%02d.img',s);
	res{s} = CM_run_mvpa_across_subjects(subj_array, trte_str, 'roiName', mask_to_use);

end
    



