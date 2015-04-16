function [h] = sl_mvpa_group_vis(analysis_name, thresh, saveSPM)
% function [h] = sl_mvpa_group_vis(analysis_name, thresh, saveSPM)
% For some searchlight mvpa retrive (or create) a group-level results 
% structure and save out niftis of average aucs threshholding 
% according to thresh

FUNC_DIR = '/hsgs/projects/awagner/jtboyd/CM_ret/Data/Functional/';
SL_DIR = fullfile(FUNC_DIR,'sl_mvpa');
analysis_dir = fullfile(SL_DIR, analysis_name);
group_file = fullfile(analysis_dir, 'group_summary.mat');

SUBJ_EG_FILE = fullfile(FUNC_DIR, 'CM001/mvpa/CM001_SEPT09_MVPA_MASK_resliced4mm.nii_s8mm_wa.mat');

% load for the corresponding group summary file. if it doesn't exist, create it
if exist(group_file, 'file')
	g = load(group_file);
else
	g.group = sl_mvpa_postProc(analysis_name,0);
end

% find voxels who showed values significantly different from chance (where significance is defined by thresh)
sig_ix = find(g.group.p<thresh);
supra_thresh_ix_pos = find(g.group.p<thresh & g.group.mean_auc > .5);
supra_thresh_ix_neg = find(g.group.p<thresh & g.group.mean_auc < .5);
%mean_auc_above_chance = nan(size(g.group.mean_auc));
%mean_auc_below_chance = nan(size(g.group.mean_auc));
mean_auc_diff_chance = nan(size(g.group.mean_auc));
%mean_auc_above_chance(supra_thresh_ix_pos) = g.group.mean_auc(supra_thresh_ix_pos);
%mean_auc_below_chance(supra_thresh_ix_neg) = g.group.mean_auc(supra_thresh_ix_neg);
mean_auc_diff_chance(sig_ix) = g.group.mean_auc(sig_ix);
% remove .5, so that we can plot negative and positive values in caret
mean_auc_diff_chance_minuspt5 = mean_auc_diff_chance - .5;
% for now, let's load a relevant subject classification struct
subj_eg = load(SUBJ_EG_FILE);
subj_eg.subj.header.id='group';
subj_eg = initset_object(subj_eg.subj, 'pattern', sprintf('auc_diff_chance_p%s_minuspt5',num2str(thresh)), mean_auc_diff_chance_minuspt5','masked_by','SEPT09_MVPA_MASK_resliced4mm.nii')
%view_pattern_overlay(subj_eg,'SEPT09_MVPA_MASK_resliced4mm.nii','above_thresh_auc')
subj_eg.patterns{end}.header.vol = subj_eg.patterns{end-1}.header.vol{1}
subj_eg.patterns{end}.header.vol.fname = fullfile(analysis_dir,sprintf('mean_auc_diff_chance_p%s_minuspt5.nii',num2str(thresh)))
write_to_spm(subj_eg,'pattern',sprintf('auc_diff_chance_p%s_minuspt5',num2str(thresh)))

