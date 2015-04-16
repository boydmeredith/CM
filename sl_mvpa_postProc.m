function [group] = sl_mvpa_postProc(analysis_name, save)
S_ARR = [1, 3:10, 12:24];
SL_DIR = '/hsgs/projects/awagner/jtboyd/CM_ret/Data/Functional/sl_mvpa';
analysis_dir = fullfile(SL_DIR, analysis_name);

for iSub = S_ARR
    this_res_file = fullfile(analysis_dir, sprintf('CM%03d.mat',iSub))
    
    % check to make sure that there is a result saved for this subject,
    % otherwise warn the user
    if exist(this_res_file, 'file')
        this_res = load(this_res_file);
        
        % if the group variables have not been created, create them
        if ~exist('group','var')
		group.eg_res = this_res.res; %keep an example results structure around in case it's useful
	    group.scram   = nan(length(S_ARR), length(this_res.res.mat));
            group.unscram = nan(length(S_ARR), length(this_res.res.mat));
        end
    else
        warning(sprintf('there is no results file saved for subject CM%03d', iSub));
        continue;
    end
    
    group.unscram(iSub, :) = nanmean(this_res.res.auc,2);
    group.scram(iSub, :) = nanmean(this_res.res.auc_scram_regs,2);
end

[group.h, group.p, group.ci, group.tstats] = ttest(group.scram, group.unscram);
[group.vspt5.h group.vspt5.p] = ttest(group.unscram, .5)
group.mean_auc = nanmean(group.unscram);
group.sd_auc = nanstd(group.unscram);
group.max_auc = max(nanmean(group.unscram));
group.min_auc = min(nanmean(group.unscram));
if save
  save(fullfile(analysis_dir, 'group_summary.mat'),'group');
end
% what should go here, you ask? something to set up a visualization of the
% average aucs across all the subjects, of course!
