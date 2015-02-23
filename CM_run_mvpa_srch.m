function [subj res] = CM_run_mvpa_srch(subjArray, classification_name, proclus, varargin)
% function [subj res] = CM_run_mvpa_srch(subjArray, classification_name, proclus)
% Runs a searchlight classification analysis on the countermeasures dataset
% Args:
%   subjArray: array of subject numbers to include in the analysis
%   classification_name: will put results in file of this name
%   proclus: allows setting up paths approriately for local machine vs proclus
%
% Returns:
%   Map of highest classification performance

% ID numbers for all subjects in this study
ALL_SUBS = [1,3:10,12:26];
% this selector will be set up by cm_create_sl_selectors and its
% name will be passed to feature_select
SL_SELECTOR_TO_USE = 'trEx_teCm_sl';

% if specific subjects are not specified to classify
if isempty(subjArray)
    subjArray=ALL_SUBS;
end

res = [];
% initialize some classifier and searchlight specific parameter variables
class_args.train_funct_name = 'train_pLR';
class_args.test_funct_name = 'test_pLR';
class_args.penalty=1;

scratch.class_args = class_args;
scratch.perfmet_funct = 'perfmet_auc';
scratch.perfmet_args = struct([]);

statmap_srch_arg.obj_funct = 'statmap_classify';
statmap_srch_arg.scratch = scratch;

% iterate through each subject performing classification
for iSub = subjArray
    % load relevant mvpa and experiment parameters
    [expt, classArgs, par] = CM_mvpa_params(iSub, 'ret', proclus);
    expt = acceptUserInput(expt,varargin);
    % create necessary strings
    expt.sl_selector_to_use = SL_SELECTOR_TO_USE;
    subjId = sprintf(expt.subjIdFormat,iSub);
    expt.saveName = cm_make_sl_savename(expt);
    expt.impMapDirStr = fullfile(expt.dir, '/mvpa_results/importance_maps/', expt.saveName);
    expt.roiFname = fullfile(expt.dir,'Masks',expt.roiName);
    thisMvpaDir = fullfile(expt.dir, subjId, expt.mvpaDirStr);
    expt.subjFname = fullfile(thisMvpaDir, [subjId '_' expt.roiName '_s8mm_wa.mat']);
    
    
    % training and testing on the same TR, whether you like it or not
   % assert(isequal(expt.trWeights_train, expt.trWeights_test));
    expt.trWeights = expt.trWeights_train;
    
    % load or create the subj structure and output a summary
    if ~exist(expt.subjFname,'file')
        subj = CM_mvpa_load_and_preprocess_raw_data(subjId, expt, nRuns, 1);
    else
        subj = load(expt.subjFname,'subj');
        subj=subj.subj;
    end
    summarize(subj);
    
    % load onsets and names of scan files
    onsetsFile = fullfile(thisMvpaDir, expt.onsetsFname);
    expt.ons = load(onsetsFile);
    expt.scanfiles = vertcat(par.swascanfiles.(par.task));
    
    % condense onsets and patterns
    fprintf('\n\nPrepping onsets and patterns for CM%03d...',iSub);
    subj = cm_condense_onsets_and_patterns(subj, expt);
    summarize(subj);
    
    active_trials = find(sum(get_mat(subj,'regressors','conds')));
    actives_selector = zeros(1,size(get_mat(subj,'regressors','conds'),2));
    actives_selector(active_trials) =1;
    subj = initset_object(subj,'selector','conditions_of_interest',actives_selector);
    
  % create selectors that will determine training and testing scheme
    fprintf('\n\nPrepping selectors for CM%03d...',iSub);
    subj = cm_create_sl_selectors(subj);
    
    summarize(subj);
    
    subj_prebalancing = subj;
    
    for iResIteration = 1:expt.num_results_iter
        
        subj = subj_prebalancing;
        new_active_trials = [];
        new_actives_selector = zeros(1,size(subj.regressors{1}.mat,2));
        
        if expt.equate_number_of_trials_in_cond_1_and_2
            subj = create_balanced_xvalid_selectors_searchlight(subj, 'conds',SL_SELECTOR_TO_USE);
        end
        
        nSelectors = length(find_group(subj,'selector',[SL_SELECTOR_TO_USE '_bal']));
        % assert only one selector, since when we have multiple selectors
        % to use it becomes difficult to figure out how to set up selectors
        % ... although, if the active_trials is just relevant for zscoring,
        % then why don't we just use the balanced selectors
        assert(nSelectors==1);
        for iSelector = 1:nSelectors
            new_active_trials = horzcat(new_active_trials, ...
                find(ismember(subj.selectors{end-nSelectors+iSelector}.mat,[1 2 3])));
        end
        new_active_trials = unique(new_active_trials);
        new_actives_selector(new_active_trials)=1;
        subj = initset_object(subj,'selector',...
            'conditions_of_interest_bal_within_runs', new_actives_selector);
        active_trials = new_active_trials;
        
        if expt.perform_second_round_of_zscoring
            display('Performing second round of z-scoring');
            pat_ix = get_number(subj,'pattern','epi_d_hp_z_condensed');
            subj.patterns{pat_ix}.mat(:,active_trials) = zscore(subj.patterns{pat_ix}.mat(:,active_trials)')';
        end
        
             
        % create spheres to use for classification. given some radius,
        % produce a matrix of nvox_in_mask x nvox_in_sphere where row i
        % contains the voxels in the sphere centered at voxel i
        subj.adj_sphere = create_adj_list(subj,expt.roiName,'radius',expt.srch_radius);
        statmap_srch_arg.adj_list = subj.adj_sphere;
        
        
        
        subj = JR_scramble_regressors(subj,'conds',[SL_SELECTOR_TO_USE '_bal_1'],'conditions_of_interest_bal_within_runs','conds_scrambled');
       
        
        % run searchlight classification
        subj = feature_select(subj, ...
            'epi_d_hp_z_condensed', ... %data
            'conds', ... % binary regs
            SL_SELECTOR_TO_USE, ... % selector
            'statmap_funct', 'statmap_searchlight', ... %function
            'statmap_arg',statmap_srch_arg, ...
            'new_map_patname', 'epi_d_hp_z_condensed_srch', ...
            'thresh', []);
        
        % rerun with the scrambled regressors
        subj = feature_select(subj, ...
            'epi_d_hp_z_condensed', ... %data
            'conds_scrambled', ... % binary regs
            SL_SELECTOR_TO_USE, ... % selector
            'statmap_funct', 'statmap_searchlight', ... %function
            'statmap_arg',statmap_srch_arg, ...
            'new_map_patname', 'epi_d_hp_z_condensed_srch_scrambled', ...
            'thresh', []);
        
               
        % visualize 1st resulting searchlight pattern
        figure;
        subj = load_spm_mask(subj,'wholebrain',fullfile(expt.dir, '/Masks/SEPT09_MVPA_MASK_resliced4mm.nii'));
        if ~proclus
		view_pattern_overlay(subj,'wholebrain','epi_d_hp_z_condensed_srch_1',...
            'over_cmap','redblue','autoslice',true)
    	end
        % populate res structure to save
        if iResIteration == 1
            res = subj.patterns{end};
            res.masks =  subj.masks;
            res.auc = nan(subj.patterns{end}.matsize(1), expt.num_results_iter);
            res.auc_scram_regs = nan(subj.patterns{end}.matsize(1), expt.num_results_iter);
            res.parameters = expt;
        end
        res.auc(:,iResIteration) = get_mat(subj,'pattern','epi_d_hp_z_condensed_srch_1');
        res.auc_scram_regs(:,iResIteration) = get_mat(subj,'pattern','epi_d_hp_z_condensed_srch_scrambled_1');
        
%         % hack together header information for the new pattern, so that we
%         % can save out a Nifti corresponding to the new info
%         assert(length(subj.patterns) == 6) % I'm being lazy in how the header gets defined, so let's make sure to assert that I change it in the future
%         subj.patterns{6}.header.vol = subj.patterns{5}.header.vol{1}
%         subj.patterns{6}.header.vol.fname = ...
%             sprintf(fullfile(expt.dir, '/CM%03d/test_searchlight_output_%s.nii'),iSub,expt.roiName);
%         write_to_spm(subj,'pattern','epi_d_hp_z_condensed_srch_1');
        
        
    end
    
    saveDir = fullfile(expt.dir, 'sl_mvpa', classification_name);
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    save(fullfile(saveDir, sprintf('CM%03d.mat',iSub)),'res');
    close all;
end

% if ~exist(expt.groupMvpaDir)
%     mkdir(expt.groupMvpaDir);
% end
% save(fullfile(expt.group_mvpa_dir, expt.saveName),'res', 'subj');


end

function name = cm_make_sl_savename(expt)
condNamesStr = strrep(strjoin(expt.condNames),',','V');
trTeStr = strrep(num2str(expt.which_traintest),'  ','_');
trWStr_tr = strrep(num2str(expt.trWeights_train),' ','_');
roiStr = regexprep(expt.roiName,'(\w*).*','$1');
name = sprintf('srchlight_conds_%s_trTe%s_trW%s_roi%s',condNamesStr,trTeStr,trWStr_tr,roiStr);

end

function subj = cm_condense_onsets_and_patterns(subj, expt)
% Condenses the loaded onsets and patterns to remove timepoints of
% non-interest
condensed_runs = [];
nExamples = [];

names = expt.ons.names;
onsets = expt.ons.onsets;

nScanFiles = length(expt.scanfiles);
nRuns = nScanFiles/expt.numTpPerRun;
nPatts = size(subj.patterns{end}.mat,2)
nCondsTotal = size(onsets,2);
assert(nScanFiles == nPatts);

% find the
expt.condCols = makeCondCols(expt, names);
nCondsToClassify = length(expt.condNames);
condOnsets = zeros(nCondsToClassify, nPatts);
i = 1;
for iCond = expt.condCols;
    % record number of examples available for each condition
    nExamples(i) = length(onsets{iCond});
    % turn the onsets into a value that corresponds with volume acquisition
    % number, rather than time passed
    time_idx = onsets{iCond}/expt.trSecs + 1;
    % put 1's in the condOnsets to create regressors
    condOnsets(i, time_idx) = 1;
    i = i + 1;
end
% visualize the trials associated with each condition
subplot(1,2,1);
imagesc(condOnsets');

% identify rest timepoints (defined here as timepoints where neither
% condition of interest is relevant
restTp = (sum(condOnsets,1) == 0);
% create condensed regressors, which exlude timepoints of noninterest
condensedCondRegs = condOnsets(:,~restTp);
subplot(1,2,2);
imagesc(condensedCondRegs');
runs = get_mat(subj,'selector','runs');
all_trials = sum(condOnsets,1);
% if we're training on the explicit data and testing on the countermeasures
% we should make sure to accomodate training and testing on different trs
if strcmp(expt.sl_selector_to_use, 'trEx_teCm_sl') 
ex_trials = runs < 5 & all_trials;
cm_trials = runs >= 5 & all_trials;
	for trNum = expt.trsToAverageOver
	    data_by_tr_train(trNum,:,:) = expt.trWeights_train(trNum)*subj.patterns{end}.mat(:,find(ex_trials)+trNum-1);
	    data_by_tr_test(trNum,:,:) = expt.trWeights_test(trNum)*subj.patterns{end}.mat(:,find(cm_trials)+trNum-1);
	    data_by_tr = cat(3,data_by_tr_train, data_by_tr_test);
	end
else
	for trNum = expt.trsToAverageOver
	    data_by_tr(trNum,:,:) = expt.trWeights(trNum)*subj.patterns{end}.mat(:,find(all_trials)+trNum-1);
	end
end
temporally_condensed_data = squeeze(sum(data_by_tr(expt.trsToAverageOver,:,:),1));
clear data_by_tr;

if expt.remove_outlier_trials
    mean_across_voxels = mean(temporally_condensed_data,1);
    z_mean_across_voxels = zscore(mean_across_voxels);
    upper_outliers = find(z_mean_across_voxels> expt.remove_outlier_trials);
    lower_outliers = find(z_mean_across_voxels< -1 * expt.remove_outlier_trials);
    all_outliers = union(upper_outliers,lower_outliers)
    condensedCondRegs(:,all_outliers) = 0;
    fprintf('removing outliers >%i standard deviations from the mean (trials %s)',expt.remove_outlier_trials,all_outliers)
end

% create condensed runs be removing rest timepoints
condensed_runs = subj.selectors{1}.mat(~restTp);

% add condensed images to subj
subj = initset_object(subj, 'regressors', 'conds', condensedCondRegs, 'condnames', expt.condNames);
subj = duplicate_object(subj, 'pattern', 'epi_d_hp_z','epi_d_hp_z_condensed');
subj = set_mat(subj, 'pattern', 'epi_d_hp_z_condensed', temporally_condensed_data,'ignore_diff_size',true);
subj = add_history(subj,'pattern','epi_d_hp_z_condensed','Pattern created created by JR custom code');
subj = remove_mat(subj,'pattern','epi_d_hp_z');
summarize(subj);

% add condensed runs variable to subj
subj.selectors{1}.mat = condensed_runs;
subj.selectors{1}.matsize = size(condensed_runs);

end


function subj = cm_create_sl_selectors(subj)
FIRST_CM_RUN = 5;
CM_RUNS = 5:8;
runs = get_mat(subj,'selector','runs');
nRuns = max(runs);
nTimepoints = length(runs);

runs_xval_sl = ones(nRuns, nTimepoints);
% % set up train explicit test countermeasures selector by setting all explicit trials to be training examples (1s), 3 of 4 CM trials are set as searchlight training runs (3s), and the remaining CM trial is set as a final generalization test trial (2)
% trEx_teCm_sl = ones(4, nTimepoints);
% cm_run_ix = find(ismember(runs,CM_RUNS));
% trEx_teCm_sl(:, cm_run_ix) = 3;
% for r = FIRST_CM_RUN:nRuns
%     r_ix = r-FIRST_CM_RUN+1;
%     cur_final_test_run = find(runs == nRuns-r_ix+1);
%     trEx_teCm_sl(r_ix, cur_final_test_run) = 2;
%     cur_name = sprintf('trEx_teCm_sl_%i',r_ix);
%     subj = initset_object(subj, 'selector', cur_name, ...
%         trEx_teCm_sl(r_ix,:),'group_name','trEx_teCm_sl');
% end
% 
% imagesc(trEx_teCm_sl);
% set(gca, 'CLim', [1 3]);
% colorbar;

trEx_teCm_sl = ones(1,nTimepoints);
cm_run_ix = find(ismember(runs,CM_RUNS));
trEx_teCm_sl(:, cm_run_ix) = 3;
trEx_teCm_sl(:,subj.selectors{end}.mat==0) = 0;
subj = initset_object(subj, 'selector', 'trEx_teCm_sl_1', ...
         trEx_teCm_sl,'group_name','trEx_teCm_sl');

% set up xval selector
for r = 1:nRuns
    cur_final_test_run = find(runs==r);
    runs_xval_sl(r, cur_final_test_run) = 2;
end

%imagesc(runs_xval_sl);
%set(gca, 'CLim', [1 3]);
%colorbar

for r = 1:nRuns
    cur_searchlight_gen_run = find(runs == nRuns-r+1);
    runs_xval_sl(r, cur_searchlight_gen_run) = 3;
    runs_xval_sl(r,subj.selectors{end}.mat==0) = 0;
    cur_name = sprintf('runs_xval_sl_%i',r);
    subj = initset_object(subj, 'selector', cur_name, ...
        runs_xval_sl(r,:), ...
        'group_name', 'runs_xval_sl' );
end

%imagesc(runs_xval_sl)
%colorbar



end




