function []= CM_run_mvpa_v3_searchlight(subj_array, condNames, radius, penalty, which_traintest, nickname)

% Workflow script for running binary MVPA classification on Jesse Rissman's
% face recognition memory dataset. Requires Princeton MVPA toolbox
% functions to be located in your path.

% example usage:  CM_run_mvpa_v3({'CM001'},'EX_hits','EX_CRs',3,10)

% Inputs
%
% subj_array: array of subject id #'s (e.g., [1 2 4 5 8]) to run analysis on; these #'s get converted to strings below (e.g., 's111', 's112', 's114', etc.)
%
% condition1: name of condition1 (e.g., 'Hits')
% condition2: name of condition2 (e.g., 'CRs')
%
% radius:  searchlight radius in voxels (e.g., 3 voxel radius --> 123 voxel
% sphere size)


%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)

flags.num_full_iter = 1; % number of times to run the entire classification process, including feature selection
flags.num_results_iter = 10; % number of times to run the post-feature selection classification process (select subset of the data and train/test classifier)
flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
flags.equate_number_of_trials_in_cond_1_and_2 = 1; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
flags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
flags.anova_nVox_thresh = 0; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
flags.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification (because trial selection and TR selection undoes initial z-scoring)
flags.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
flags.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
flags.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
flags.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
flags.minimum_number_of_trials = 15;
flags.generate_importance_maps = 1; % 1=generate importance maps based on classification weights (scaled by the mean univariate activity of each condition)
flags.write_data_log_to_text_file=1; % save a .txt file that logs critical classification performance data
flags.save_data_log_as_mat_file =0; % save a .mat file that logs critical classification performance data
flags.lambda = penalty; % penalty parameter for logistic regression classifier
flags.optimize_penalty_param = 0; % 1 = peformed nested cross-validation for empirical penalty optimization
flags.searchlight_radius = radius;

% specify classifier
class_args.train_funct_name = 'train_pLR';
class_args.test_funct_name = 'test_pLR';
class_args.penalty = flags.lambda;

for subNum=1:length(subj_array)  % run loop for each subject in subj array (subject ID #'s specified in short format)
    tic % record the start time
    subj_id = subj_array{subNum};
    expt_name = 'COUNTERMEASURES';
    expt_dir = '/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional/'; % top-level expt. directory where individual subject folders live
    mvpa_dir = [expt_dir '/' subj_id '/mvpa']; % mvpa directory within each subject's folder (will be created below if it doesn't exist)
    % change mask back to SEPT09_MVPA_MASK ?????
    roi_name = 'SEPT09_MVPA_MASK_resliced4mm';  % name of mask to be used for voxel selection (can be a small ROI, a whole-brain mask, or anywhere in between)
    roi_file = [expt_dir '/Masks/' roi_name '.nii']; % specific path to mask file
    data_imgs_to_use = 'raw_filenames.mat'; % .mat file containing names of all functional images (must exist for each subject; can be created by running cellstr(SPM.xY.P) on subject's SPM.mat file)
    num_TP_per_run = 256; % number of TRs per scanning run (coded for fixed value; adjust code structure if variable TR counts across runs)
    load([mvpa_dir '/onsets_mvpa.mat']); % load in your SPM-formatted onsets file (specifies onsets of each event time in seconds)
    load([mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array called raw_filenames into memory
    num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
    % specify previously saved mvpa workspace to bypass time-consuming data extraction and preprocessing (if not yet created, this will specify it's name and location)
    thisSubjFname = [mvpa_dir '/' subj_id '_' roi_name '_s8mm_wa.mat'];
    load([expt_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
    %turn the condNames array into a single string to describe the comparison
    condStr = [repmat(sprintf(['%s', '_vs_'], condNames{1:end-1}),1, ~isscalar(condNames)), sprintf('%s', condNames{end})];
    
    % Old Code expected TRs to be specified as an argument. However,
    % searchlight simply uses all TRs and weights them appropriately as follows:
    TRs_to_average_over = [1 2 3 4 5]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
    TR_weights = [0 0 .33 .34 .33];
    
    traintest =  {'EXP_only','CM_only','EXP>CM','EXP>CM_first_run','EXP>CM_last_run','CM>EXP','All_runs','EXP>CM_second_run','EXP>CM_third_run', 'All_runs_TEST_CM', 'All_runs_TEST_EX'};
    
    % specify directory names, create new directories if necessary
    dir_str = [traintest{which_traintest} '_' condStr '_pen' num2str(penalty) '_rad' num2str(radius) '_' nickname];
    searchlight_dir = [expt_dir '/mvpa_results/within_subj_class/searchlight_maps/whole_brain_searchlight_maps/' dir_str];
    
    if ~exist(searchlight_dir,'dir')
        mkdir(searchlight_dir);
    end
    
    
    if ~exist(thisSubjFname,'file') %if a workspace .mat file has not yet been created for this subject
        
        % call script to run time-consuming data preprocessing routines
        % (load pattern, detrend, high-pass filter, and z-score)
        [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, expt_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);
        save(thisSubjFname, 'subj');
    else
        load(thisSubjFname)  %load workspace
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract info about conditions from onsets file (uses SPM-based coding
    % of the onsets (in seconds) of trials in each condition of interest)
    num_conds = size(onsets,2);
    
    all_onsets = zeros(num_conds,size(subj.patterns{end}.mat,2)); % initialize all_onsets matrix as conditions x timepoints
    
    for cond = 1: num_conds
        for trial = 1: length(onsets{cond})
            time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_onsets(cond,time_idx) = 1;
        end
    end
    
    % condense all_onsets matrix from having one column per TR to having
    % once column per trial (remove TRs that have 0's for all
    % conditions, i.e. rest timepoints)
    
    condensed_regs_all = [];
    condensed_runs = [];
    trial_counter = 1;
    for i = 1: size(all_onsets,2)
        if ~isempty(find(all_onsets(:,i))) % if not a rest timepoint
            condensed_regs_all(:,trial_counter) = all_onsets(:,i);
            condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
            trial_counter = trial_counter + 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % select TRs of interest (to correspond with peak post-stim BOLD response)
    all_trials = sum(all_onsets,1); % vector of all trials
    data_by_TR(1,:,:) = TR_weights(1)*subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
    data_by_TR(2,:,:) = TR_weights(2)*subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
    data_by_TR(3,:,:) = TR_weights(3)*subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
    data_by_TR(4,:,:) = TR_weights(4)*subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
    data_by_TR(5,:,:) = TR_weights(5)*subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)
    temporally_condensed_data = squeeze(sum(data_by_TR(TRs_to_average_over,:,:),1));
    clear data_by_TR; % clean up matlab workspace to save memory
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Exclude trials determined to be outliers by custom ArtDetect script
    % Guide to outlier file cell arrays...
    % Movement thresholds: .2 .25 .3 .35 .4 .4 .5
    % Global signal thresholds: 2 2.5 3 3.5 4 4.5 5
    
    if flags.remove_artdetect_outliers == 1
        load([expt_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices
        
        m_outliers = movement_outlier_trials{flags.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
        gs_outliers = global_signal_outlier_trials{flags.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
        combined_outliers = union(m_outliers,gs_outliers);
        
        condensed_regs_all(:,combined_outliers) = 0;
        
        display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
        display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
        display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % on-the-fly outlier detection/removal;
    if flags.remove_outlier_trials ~= 0
        mean_across_voxels = mean(temporally_condensed_data,1);
        z_mean_across_voxels = zscore(mean_across_voxels);
        upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
        lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
        all_outliers = union(upper_outliers,lower_outliers)
        condensed_regs_all(:,all_outliers) = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    subj_original = subj; % backup the subj structure before entering k-loop
    condensed_regs_all_original = condensed_regs_all; % backup condensed_regs_all before entering k-loop
    x = 0;  %initialize the counter x (gets incremented during each classification run-through)
    for k = 1:flags.num_full_iter
        
        subj = subj_original;
        condensed_regs_all = condensed_regs_all_original;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define conditions of interest (EXPERIMENT-SPECIFIC CODE NEEDED TO
        % INDEX TRIALS FROM EACH CONDITION OF INTEREST)
        
        
        % define mnemonic states based on individual conditions or combinations of conditions
        
        %CONDITION NUMBERS:
        %{1='AA_Explicit_Hit',2='AA_Explicit_Miss',3='EA_Explicit_Hit',4='EA_Explicit_Miss',5='AA_Explicit_CR',6='AA_Explicit_FA',...
        %...7='EA_Explicit_CR',8='EA_Explicit_FA',9='AA_CM_Hit',10='AA_CM_Miss',11='EA_CM_Hit',12='EA_CM_Miss',...
        %...13='AA_CM_CR',14='AA_CM_FA',15='EA_CM_CR',16='EA_CM_FA',17='no_response',18='leadout_instructions';}
        
        EX_OLD = sum(condensed_regs_all(1:4,:)); %"OLD" = objective old status
        EX_NEW = sum(condensed_regs_all(5:8,:));
        EX_old = sum(condensed_regs_all([1 3 6 8],:)); %"old" = subjective old status
        EX_new = sum(condensed_regs_all([2 4 5 7],:));
        
        EX_hits = sum(condensed_regs_all([1 3],:));
        EX_misses = sum(condensed_regs_all([2 4],:));
        EX_FAs = sum(condensed_regs_all([6 8],:));
        EX_CRs = sum(condensed_regs_all([5 7],:));
        
        
        CM_OLD = sum(condensed_regs_all(9:12,:));
        CM_NEW = sum(condensed_regs_all(13:16,:));
        CM_old = sum(condensed_regs_all([9 11 14 16],:));
        CM_new = sum(condensed_regs_all([10 12 13 15],:));
        
        CM_hits = sum(condensed_regs_all([9 11],:));
        CM_misses = sum(condensed_regs_all([10 12],:));
        CM_FAs = sum(condensed_regs_all([14 16],:));
        CM_CRs = sum(condensed_regs_all([13 15],:));
        
        % Race within Tasks
        AA_EX_OLD = sum(condensed_regs_all([1 2],:));
        AA_EX_NEW = sum(condensed_regs_all([5 6],:));
        AA_EX_old = sum(condensed_regs_all([1 6],:));
        AA_EX_new = sum(condensed_regs_all([2 5],:));
        
        AA_EX_OLDandNEW = sum(condensed_regs_all([1 2 5 6],:));
        EA_EX_OLDandNEW = sum(condensed_regs_all([3 4 7 8],:));
        
        AA_EX_hits = condensed_regs_all(1,:);
        AA_EX_misses = condensed_regs_all(2,:);
        AA_EX_FAs = condensed_regs_all(6,:);
        AA_EX_CRs = condensed_regs_all(5,:);
        
        
        EA_EX_OLD = sum(condensed_regs_all([3 4],:));
        EA_EX_NEW = sum(condensed_regs_all([7 8],:));
        EA_EX_old = sum(condensed_regs_all([3 8],:));
        EA_EX_new = sum(condensed_regs_all([4 7],:));
        
        EA_EX_hits = condensed_regs_all(3,:);
        EA_EX_misses = condensed_regs_all(4,:);
        EA_EX_FAs = condensed_regs_all(8,:);
        EA_EX_CRs = condensed_regs_all(7,:);
        
        AA_CM_hits = condensed_regs_all(9,:);
        AA_CM_CRs = condensed_regs_all(13,:);
        EA_CM_hits = condensed_regs_all(11,:);
        EA_CM_CRs = condensed_regs_all(15,:);
        
        AA_CM_OLDandNEW = sum(condensed_regs_all([9 10 13 14],:));
        EA_CM_OLDandNEW = sum(condensed_regs_all([11 12 15 16],:));
        
        
        %
        % %         % Race across Tasks
        AA_OLD = sum(condensed_regs_all([1 2 9 10],:));
        AA_NEW = sum(condensed_regs_all([5 6 13 14],:));
        % %         AA_old = sum(condensed_regs_all([1 6],:));
        % %         AA_new = sum(condensed_regs_all([2 5],:));
        % %
        AA_hits = sum(condensed_regs_all([1 9],:));
        % %         AA_misses = sum(condensed_regs_all(2,:));
        % %         AA_FAs = sum(condensed_regs_all(6,:));
        AA_CRs = sum(condensed_regs_all([5 13],:));
        % %
        % %
        EA_OLD = sum(condensed_regs_all([3 4 11 12],:));
        EA_NEW = sum(condensed_regs_all([7 8 15 16],:));
        % %         EA_old = sum(condensed_regs_all([3 8],:));
        % %         EA_new = sum(condensed_regs_all([4 7],:));
        % %
        EA_hits = sum(condensed_regs_all([3 11],:));
        % %         EA_misses = sum(condensed_regs_all(4,:));
        % %         EA_FAs = sum(condensed_regs_all(8,:));
        EA_CRs = sum(condensed_regs_all([7 15],:));
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% Specify which training/testing scheme you wish you to use %%%%%%%%%%%%%%%%%
        
        
        
        switch which_traintest
            case 1  %% Train and Test on EXPLICIT runs only
                condensed_runs(condensed_runs>4)=4;  % give all countermeasures runs a run label of 4; since there will be no active_trials from these runs, cross-validation will be 4-fold (i.e., bc we specify the conditions to discriminate between (e.g., EX_hits vs. EX_CRs) and they are labeled specificall for EXP or CM, don't need to ignore the CM runs here bc no CM conds specified as active)
            case 2  %% Train and Test on COUNTERMEASURES runs only
                condensed_runs(condensed_runs<=4)=5;  % give all explicit memory runs a run label of 5; since there will be no active_trials from these runs, cross-validation will be 4-fold
                condensed_runs = condensed_runs - 4;  % subtract 4 from each run label, turning run5 into run1, run6 into run2, etc.
            case 3  %% Train on EXPLICIT runs and Test on COUNTERMEASURES runs
                condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
                condensed_runs(condensed_runs>4)=1; % data from runs 5-8 will be the first testing set
            case 4  %% Train on EXPLICIT runs and Test on FIRST COUNTERMEASURES run only
                condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
                condensed_runs(condensed_runs==5)=1; % data from run 5 will be the first testing set
                condensed_runs(condensed_runs>5)=0; % data from other CM runs ignored
            case 5  %% Train on EXPLICIT runs and Test on LAST COUNTERMEASURES run only
                condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
                condensed_runs(condensed_runs==8)=1; % data from run 8 will be the first testing set
                condensed_runs(intersect(find(condensed_runs>4), find(condensed_runs<8)))=0; % data from other CM runs ignored
            case 6  %% Train on COUNTERMEASURES runs and Test on EXPLICIT runs
                condensed_runs(condensed_runs<=4)=1; % data from runs 1-4 will be the first testing set
                condensed_runs(condensed_runs>4)=2; % data from runs 5-8 will be the first training set
            case 7  %% Train and Test on ALL RUNS (across EXP and CM)
                condensed_runs = 1+mod(condensed_runs-1,4);  % training runs include a EXP and a CM run, test on remaining EXP and CM runs; cross-validation will therefore be 4-fold
            case 8  %% Train on EXPLICIT runs and Test on SECOND COUNTERMEASURES run only
                condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
                condensed_runs(condensed_runs==6)=1; % data from run 6 will be the first testing set
                condensed_runs(condensed_runs==5)=0; % data from other CM runs ignored
                condensed_runs(condensed_runs>6)=0; % data from other CM runs ignored
            case 9  %% Train on EXPLICIT runs and Test on THIRD COUNTERMEASURES run only
                condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
                condensed_runs(condensed_runs==7)=1; % data from run 7 will be the first testing set
                condensed_runs(intersect(find(condensed_runs>4), find(condensed_runs<7)))=0; % data from other CM runs ignored
                condensed_runs(condensed_runs==8)=0; % data from other CM runs ignored
            case 10 %%Train on first two explicts runs, and first two CM runs, test on second two CM runs
                condensed_runs(mod(condensed_runs-1, 4)<2) = 2;
                condensed_runs(condensed_runs>6) = 1;
                condensed_runs(intersect(find(condensed_runs>2), find(condensed_runs<5)))=0;
                
            case 11 %%Train on first two explicts runs, and first two CM runs, test on second two explicit runs runs
                condensed_runs(mod(condensed_runs-1, 4)<2) = 2;
                condensed_runs(condensed_runs>6) = 0;
                condensed_runs(intersect(find(condensed_runs>2), find(condensed_runs<5)))=1;
        end
        num_runs = length(unique(nonzeros(condensed_runs)));
        
        if which_traintest==3 || which_traintest==4 || which_traintest==5 || which_traintest==8 || which_traintest==9
            % code the CM_OLD trials as examples of condition 1 (e.g., Hits) and the CM_NEW trials as examples of condition 2 (e.g., CRs)
            eval([condNames{1} ' = ' condNames{1} ' + CM_hits;']);
            eval([condNames{2} ' = ' condNames{2} ' + CM_CRs;']);
        elseif which_traintest==6
            eval([condNames{1} ' = ' condNames{1} ' + EX_hits;']);
            eval([condNames{2} ' = ' condNames{2} ' + EX_CRs;']);
        end
        
        %%this line is a gross hack!!!%%
        if which_traintest >=3 && which_traintest ~= 7
            num_testing_sets = 1;  %% only summarize and record data from COUNTERMEASURES runs (1st cross-validation testing set)
        else
            num_testing_sets = num_runs;
        end
        
        % NOTE: this variable is specified differently in Jesse's
        % searchlight script
        num_runs = length(unique(nonzeros(condensed_runs)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %assign conditions to train/test classifier on
        condensed_regs_of_interest = [];
        eval(['condensed_regs_of_interest(1,:) = ' condNames{1} ';'])
        eval(['condensed_regs_of_interest(2,:) = ' condNames{2} ';'])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set up regressors and selectors
        
        % initialize regressors object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
        subj = set_objfield(subj,'regressors','conds','condnames',condNames);
        
        % add new condensed activation pattern
        subj = duplicate_object(subj,'pattern','epi_d_hp_z','epi_d_hp_z_condensed');
        subj = set_mat(subj,'pattern','epi_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);
        zhist = sprintf('Pattern ''%s'' created by JR custom code','epi_d_hp_z_condensed');
        subj = add_history(subj,'pattern','epi_d_hp_z_condensed',zhist,true);
        
        % clean up workspace to save RAM
        subj = remove_mat(subj,'pattern','epi_d_hp_z'); % remove original uncondensed pattern (full timeseries)
        
        % update run vector to condensed format
        subj.selectors{1}.mat = condensed_runs;
        subj.selectors{1}.matsize = size(condensed_runs);
        
        % "activate" only those trials of interest (from regs_of_interest) before creating cross-validation indices
        active_trials = find(sum(condensed_regs_of_interest)); % find those trials where the sum (across columns) is not zeros
        
        actives_selector = zeros(1,size(condensed_regs_all,2)); % create new selector object; begin by intializing vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        
        subj = init_object(subj,'selector','conditions_of_interest'); %initialize new selector object called 'conditions_of_interest'
        subj = set_mat(subj,'selector','conditions_of_interest',actives_selector); % set this selector object to have the values specified by actives_selector
        
        
        subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest', 'ignore_jumbled_runs', true,'ignore_runs_zeros',true); % create cross-validation indices (new selector group), using only the the trials specified by actives_selector
        
        
        
        % run feature selection ANOVA: specify pvalue (if desired)
        statmap_arg.use_mvpa_ver = true;
        if flags.anova_p_thresh ~= 1
            
            subj = JR_feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        % run feature selection ANOVA: specify #of voxels (if desired)
        if flags.anova_nVox_thresh ~=0
            
            subj = JR_feature_select_top_N_vox(subj,'epi_d_hp_z_condensed','conds','runs_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subj_prebalancing = subj; % backup the subj structure before entering n-loop
        active_trials_prebalancing = active_trials; % backup condensed_regs_all before entering n-loop
        
        for resultsIterNum = 1: flags.num_results_iter
            
            % restore original versions before continuing with next loop iteration
            active_trials = active_trials_prebalancing;
            subj = subj_prebalancing;
            
            % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %   NOTE: THIS PIECE OF SCRIPT WAS TO TRAIN/TEST WITHIN RUNS, SO
            % REPLACING WITH OLD LANGUAGE
            % %     % balance the number of Class A and Class B trials within each run
            % %             % (prevents classifier from becoming biased to predict one
            % %             % condition much more frequently than the other)
            % %
            % %             subj = create_balanced_xvalid_selectors(subj,'conds','runs_xval'); % creates a new 'runs_xval_bal' selector group
            % %
            % %             % after running the create_balanced_xvalid_selectors function,
            % %             % which excludes random subsets of trials from each run to
            % %             % balance the classes, we need to get a new list of the 'active trials'
            % %             new_active_trials =[];
            % %             for rr = 1:num_runs
            % %                 new_active_trials = horzcat(new_active_trials, find(subj.selectors{end-num_runs+rr}.mat==2));
            % %             end
            % %             new_actives_selector= zeros(1,size(condensed_regs_all,2));
            % %             new_actives_selector(new_active_trials)=1;
            % %             subj = init_object(subj,'selector','conditions_of_interest_bal_within_runs'); %initialize selector object
            % %             subj = set_mat(subj,'selector','conditions_of_interest_bal_within_runs',new_actives_selector);
            % %             active_trials = new_active_trials;
            % %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % %  OLD LANGUAGE STARTS HERE...
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % balance the number of Class A and Class B trials within each run
            % (prevents classifier from becoming biased to predict one
            % condition much more frequently than the other)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update the number of trials in each condition and record the
            % new set of indices
            cond1_trials = find(condensed_regs_of_interest(1,active_trials));
            cond2_trials = find(condensed_regs_of_interest(2,active_trials));
            num_cond1 = length(cond1_trials);
            num_cond2 = length(cond2_trials);
            
            % break and move on to next subject if either condition has
            % an insufficient number of trials
            min_trials = min(num_cond1,num_cond2);
            if min_trials<flags.minimum_number_of_trials
                display(['Only ' num2str(min_trials) ' trials available for ' subj_id '. Moving on to next subject...'])
                return;
            end
            
            if num_cond1 > num_cond2
                rand_array = rand(1,num_cond1);
                [sorted inds]= sort(rand_array);
                trials_to_cut = cond1_trials(inds(1:num_cond1-num_cond2));
                %condensed_regs_of_interest(1,active_trials(trials_to_cut)) = 0;
                active_trials(trials_to_cut) = [];
                
                display([num2str(length(trials_to_cut)) ' trials cut from ' condNames{1}]);
            elseif num_cond1 < num_cond2
                rand_array = rand(1,num_cond2);
                [sorted inds]= sort(rand_array);
                trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
                %condensed_regs_of_interest(2,active_trials(trials_to_cut)) = 0;
                active_trials(trials_to_cut) = [];
                display([num2str(length(trials_to_cut)) ' trials cut from ' condNames{2}]);
            else
                display('Trial numbers are already balanced');
            end
        end
        
        actives_selector = zeros(1,size(condensed_regs_all,2)); % intialize vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % balance the number of Class A and Class B trials within each run
        % (prevents classifier from becoming biased to predict one
        % condition much more frequently than the other)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update the number of trials in each condition and record the
        % new set of indices
        
        
        cond1_trials = intersect(find(condensed_regs_of_interest(1,:)),active_trials);
        cond2_trials = intersect(find(condensed_regs_of_interest(2,:)),active_trials);
        
        num_cond1 = length(cond1_trials);
        num_cond2 = length(cond2_trials);
        
        display([num2str(num_cond1) ' trials in condition ' condNames{1}])
        display([num2str(num_cond2) ' trials in condition ' condNames{2}])
        
        num_xvalid_bins = 10;
        num_runs = num_xvalid_bins;
        
        run_labels = repmat([1:num_xvalid_bins],1,ceil(num_cond1/num_xvalid_bins));
        run_labels = run_labels(1:num_cond1); %truncate the length of this vector
        
        shuffled_run_inds_for_cond1_trials = shuffle(run_labels);
        shuffled_run_inds_for_cond2_trials = shuffle(run_labels);
        
        
        subj = duplicate_object(subj,'selector','runs','runs_final');
        subj.selectors{end}.mat(cond1_trials)=shuffled_run_inds_for_cond1_trials;
        subj.selectors{end}.mat(cond2_trials)=shuffled_run_inds_for_cond2_trials;
        
        subj.selectors{end}.mat(subj.selectors{end}.mat>num_xvalid_bins) = 1;  %replace any extra xvalid labels with 1
        
        subj = init_object(subj,'selector','conditions_of_interest_final'); %initialize selector object
        subj = set_mat(subj,'selector','conditions_of_interest_final',actives_selector);
        
        subj = create_xvalid_indices(subj,'runs_final','actives_selname','conditions_of_interest_final','ignore_jumbled_runs','true','ignore_runs_zeros',true);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if flags.perform_second_round_of_zscoring == 1  % z-score the data prior to classification (active trials only)
            subj.patterns{5}.mat(:,active_trials) = zscore(subj.patterns{5}.mat(:,active_trials)')';
            display('Performing second round of z-scoring')
        end
        
        %% begin searchlight-specific code %%
        num_selectors = length(subj.selectors);
        for q = num_selectors-num_runs+1:num_selectors;
            subj.selectors{q}.mat(subj.selectors{q}.mat==2)=3; % change 2's to 3's
        end
        
        subj.patterns{end}.mat = double(subj.patterns{end}.mat);
        
        
        searchlight_radius = flags.searchlight_radius;
        subj.adj_sphere = create_adj_list(subj,roi_name,'radius',searchlight_radius);
        
        %specify training/testing functions for within-sphere classification
        class_args.train_funct_name = 'train_GNB';
        class_args.test_funct_name = 'test_GNB';
        %class_args.penalty = 1;
        
        scratch.class_args = class_args;
        scratch.perfmet_funct = 'perfmet_maxclass';
        scratch.perfmet_args = struct([]);
        
        statmap_srch_arg.adj_list = subj.adj_sphere;
        statmap_srch_arg.obj_funct = 'JR_statmap_classify_GNB_inlined';
        %statmap_srch_arg.obj_funct = 'JR_statmap_classify';
        statmap_srch_arg.scratch = scratch;
        
        subj = JR_feature_select_searchlight( ...
            subj, ...
            'epi_d_hp_z_condensed', ... % data
            'conds', ... % binary regs (for GNB)
            'runs_final_xval', ... % selector
            'statmap_funct','JR_statmap_searchlight', ... % function
            'statmap_arg',statmap_srch_arg, ...
            'new_map_patname','epi_d_hp_z_condensed_srch', ...
            'thresh',[]);
        
        % create masks from the statmaps, by picking the best N values (e.g., 500) in each
        %statmap to use for classification on the remaining run of test trials
        %NOTE: keeps top N voxels; not top N spheres
        %         subj = create_sorted_mask( ...
        %             subj,'spiral_d_hp_z_condensed_srch', ...
        %             'spiral_d_hp_z_condensed_srch_1000',1000, ...
        %             'descending',true);
        
        
        %% WRITE OUT MEAN SEARCHLIGHT MAP TO .IMG FILE
        
        % initialize empty matrices for writing searlight results
        sl_top100_pct(:,resultsIterNum) = zeros(length(subj.acts_vector),1);
        sl_top75_pct(:,resultsIterNum) = zeros(length(subj.acts_vector),1);
        sl_top50_pct(:,resultsIterNum) = zeros(length(subj.acts_vector),1);
        sl_top25_pct(:,resultsIterNum) = zeros(length(subj.acts_vector),1);
        sl_AUC(:,resultsIterNum) = zeros(length(subj.acts_vector),1);
        desireds_vector = [];
        
        % compute desireds vector
        for nn = 1:num_runs
            desireds_vector = [desireds_vector subj.regressors{1}.mat(1,subj.selectors{end-num_runs+nn}.mat==3)];
        end
        desireds_vector(desireds_vector==0)=2;
        
        %initialize hit_rate and fa_rate
        hit_rate = zeros(length(desireds_vector),1);
        fa_rate = zeros(length(desireds_vector),1);
        
        % sort by absolute value of classifier "confidence"
        for s = 1:length(subj.acts_vector)
            if ~isnan(subj.acts_vector{s}(1)) % some voxels have all NaNs for their acts_vector, so the corrects_vector for these voxels is meaningless
                [abs_sorted_diffs ind] = sort(abs(subj.acts_vector{s}),2,'descend');
                abs_correct_sorted = subj.corrects_vector{s}(ind);
                num_trials = length(abs_correct_sorted);
                
                sl_top100_pct(s,resultsIterNum)=mean(abs_correct_sorted(1:ceil(num_trials*1)));
                sl_top75_pct(s,resultsIterNum)=mean(abs_correct_sorted(1:ceil(num_trials*.75)));
                sl_top50_pct(s,resultsIterNum)=mean(abs_correct_sorted(1:ceil(num_trials*.50)));
                sl_top25_pct(s,resultsIterNum)=mean(abs_correct_sorted(1:ceil(num_trials*.25)));
                
                
                %%%%  Compute AUC of sphere %%%%
                
                
                [sorted_diffs ind] = sort(subj.acts_vector{s},2,'descend');
                correct_sorted = subj.corrects_vector{s}(ind);
                desireds_sorted = desireds_vector(ind);
                
                class_A_inds = find(desireds_sorted == 1);
                class_B_inds = find(desireds_sorted == 2);
                num_class_A = length(class_A_inds);
                num_class_B = length(class_B_inds);
                
                
                % create continuous ROC function
                for i = 1:length(sorted_diffs);
                    hit_rate(i) = count(class_A_inds<=i) / num_class_A;
                    fa_rate(i) = count(class_B_inds<=i) / num_class_B;
                end
                
                auc_overall = auroc(hit_rate,fa_rate);
                sl_AUC(s,resultsIterNum)=auc_overall;
            end
        end
    end
end

%average across results interations and multiply by 100 to convert from
%proportion to percentage
sl_top100_pct_mean = mean(sl_top100_pct,2)*100;
sl_top75_pct_mean = mean(sl_top75_pct,2)*100;
sl_top50_pct_mean = mean(sl_top50_pct,2)*100;
sl_top25_pct_mean = mean(sl_top25_pct,2)*100;
sl_AUC_mean = mean(sl_AUC,2);

critical_accuracy = binoinv(.95,num_trials,.5)*100/num_trials;
sl_p05_binomial_thresh = zeros(size(sl_AUC_mean));
sl_p05_binomial_thresh(sl_top100_pct_mean>=critical_accuracy) = 1;

sl_map = zeros(vol_info.dim);
included_voxels = find(subj.masks{1}.mat);

vol_info.fname = [searchlight_dir '/' subj_id '_' condNames{1} '_vs_' condNames{2} '_' num2str(searchlight_radius) 'vox_radius_top_100pct_' num2str(num_trials) '_trials.img'];
sl_map(included_voxels) = sl_top100_pct_mean;
spm_write_vol(vol_info, sl_map);

vol_info.fname = [searchlight_dir '/' subj_id '_' condNames{1} '_vs_' condNames{2} '_' num2str(searchlight_radius) 'vox_radius_top_75pct_' num2str(num_trials) '_trials.img'];
sl_map(included_voxels) = sl_top75_pct_mean;
spm_write_vol(vol_info, sl_map);

vol_info.fname = [searchlight_dir '/' subj_id '_' condNames{1} '_vs_' condNames{2} '_' num2str(searchlight_radius) 'vox_radius_top_50pct_' num2str(num_trials) '_trials.img'];
sl_map(included_voxels) = sl_top50_pct_mean;
spm_write_vol(vol_info, sl_map);

vol_info.fname = [searchlight_dir '/' subj_id '_' condNames{1} '_vs_' condNames{2} '_' num2str(searchlight_radius) 'vox_radius_top_25pct_' num2str(num_trials) '_trials.img'];
sl_map(included_voxels) = sl_top25_pct_mean;
spm_write_vol(vol_info, sl_map);

vol_info.fname = [searchlight_dir '/' subj_id '_' condNames{1} '_vs_' condNames{2} '_' num2str(searchlight_radius) 'vox_radius_AUC_' num2str(num_trials) '_trials.img'];
sl_map(included_voxels) = sl_AUC_mean;
spm_write_vol(vol_info, sl_map);

vol_info.fname = [searchlight_dir '/' subj_id '_' condNames{1} '_vs_' condNames{2} '_' num2str(searchlight_radius) 'vox_radius_sig_voxels_mask' num2str(num_trials) '_trials.img'];
sl_map(included_voxels) = sl_p05_binomial_thresh;
spm_write_vol(vol_info, sl_map);



time2finish = toc/60;
display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);

% %     keep subj_array condition1 condition2 balance_per_subj_bin
keep subj_array condNames which_traintest penalty radius nickname flags class_args
end




