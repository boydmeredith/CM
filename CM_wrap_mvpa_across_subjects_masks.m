function [res] = CM_wrap_mvpa_across_subjects_masks(trte,subj_array)
%CM_wrap_mvpa_across_subjects_masks() is intended to run classification
%across subjects using pre-specified masks

if isempty(subj_array)
	subj_array = [1 3:10 12:26];
end

assert(ismember(trte,{'EXminus','CMminus','minus'}));

for s = subj_array
    	
	trte_str = sprintf('%s%02d',trte,s);
	mask_to_use = sprintf('exhitsgrcrs_AND_cmhitsgrcrs_minus%02d.img',s);
	res{s} = CM_run_mvpa_across_subjects(subj_array, trte_str, 'roiName', mask_to_use);

end
    



