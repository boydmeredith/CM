function [res] = CM_wrap_mvpa_across_subjects_masks(trte,subj_array,te_subj_array, varargin)
%function [res] = CM_wrap_mvpa_across_subjects_masks(trte,subj_array,te_subj_array, varargin) is intended to run classification
%across subjects using pre-specified masks

if isempty(subj_array)
	subj_array = [1 3:10 12:26];
end
if isempty(te_subj_array)
	te_subj_array = subj_array;
end

assert(ismember(trte,{'EXminus','CMminus','minus'}));

for s = te_subj_array
    	
	trte_str = sprintf('%s%02d',trte,s);
	mask_to_use = sprintf('exhitsgrcrs_AND_cmhitsgrcrs_minus%02d.img',s);
	res{s} = CM_run_mvpa_across_subjects(subj_array, trte_str, 'roiName', mask_to_use, 'which_traintest',1:8, varargin{:});
end
    



