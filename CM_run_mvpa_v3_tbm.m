function []= CM_run_mvpa_v3_tbm(flags, expt, class_args, subj_array, condNames, which_traintest, trnSubCondsToBal, tstSubCondsToBal, nickname)

%% specify paramters for classifier
dfltFlags.num_full_iter = 1; % number of times to run the entire classification process, including feature selection
dfltFlags.num_results_iter = 10; % number of times to run the post-feature selection classification process (select subset of the data and train/test classifier)
dfltFlags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data- can run this loop multiple times if using a stochastic (i.e. non-deterministic) classification algorithm like backpropagation neural nets
dfltFlags.equate_number_of_trials_in_cond_1_and_2 = 1; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
dfltFlags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
dfltFlags.anova_expt.nVox_thresh = 0; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
dfltFlags.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification (because trial selection and TR selection undoes initial z-scoring)
dfltFlags.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
dfltFlags.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
dfltFlags.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
dfltFlags.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
dfltFlags.generate_importance_maps = 1; % 1=generate importance maps based on classification weights (scaled by the mean univariate activity of each condition)
dfltFlags.write_data_log_to_text_file=0; % save a .txt file that logs critical classification performance data
dfltFlags.save_data_log_as_mat_file =1; % save a .mat file that logs critical classification performance data
dfltFlags.optimize_penalty_param = 0; % 1 = peformed nested cross-validation for empirical penalty optimization
dfltFlags.scramble = 0;
flags = mergestructs(flags, dfltFlags);
if isempty(class_args)
    class_args.train_funct_name = 'train_pLR';
    class_args.test_funct_name = 'test_pLR';
    class_args.penalty = 10;
end

%% do some string handling for directory organization
lbl.condStr = [repmat(sprintf(['%s', '_vs_'], condNames{1:end-1}),1, ~isscalar(condNames)), sprintf('%s', condNames{end})];
lbl.trStr = ['TRs' strrep(num2str(expt.trWeights),' ','')];
lbl.penStr = '';
if isfield(class_args,'penalty')
    if flags.optimize_penalty_param == 1
        lbl.penStr = 'OPTIMAL_pen_';
    else
        lbl.penStr = ['_pen' num2str(class_args.penalty) '_'];
    end
end
if flags.scramble ==1
    nickname = [nickname '_SCRAMBLED'];
end
lbl.dirStr = [expt.traintest{which_traintest} '_' lbl.condStr '_' class_args.train_funct_name(6:end) lbl.penStr num2str(flags.anova_expt.nVox_thresh) 'vox_' lbl.trStr '_' nickname];
lbl.impMapDirStr=[expt.dir '/mvpa_results/within_subj_class/importance_maps/' lbl.dirStr];
lbl.dataLogTxtDir=[expt.dir '/mvpa_results/within_subj_class/class_perf/' lbl.dirStr];
lbl.dataLogMatDir=[expt.dir  '/mvpa_results/within_subj_class/data_logs/' lbl.dirStr];
if flags.generate_importance_maps == 1
    if ~exist(lbl.impMapDirStr,'dir'),mkdir(lbl.impMapDirStr);end
end
if flags.write_data_log_to_text_file == 1
    if ~exist(lbl.dataLogTxtDir, 'dir'),mkdir(lbl.dataLogTxtDir);end
end
if flags.save_data_log_as_mat_file == 1
    if ~exist(lbl.dataLogMatDir, 'dir'),mkdir(lbl.dataLogMatDir);end
end

%% for every subject: preproc, classify & spit out results as specified
for subNum=1:length(subj_array)
    tic % record the start time
    %% load or create subj and temporally condense. load onsets to use as regressors and condense so that we have one per trial (instead of 1 per TR)
    subjId = sprintf(expt.subjIdFormat,subj_array(subNum));
    thisMvpaDir = fullfile(expt.dir, subjId, expt.mvpaDirStr); % mvpa directory within each subject's folder (will be created below if it doesn't exist)
    onsetsFile = fullfile(thisMvpaDir, expt.onsetsFname);
    dataImgsFile = fullfile(thisMvpaDir, expt.dataImgsToUse);
    load(onsetsFile); % load in your SPM-formatted onsets file (specifies onsets of each event time in seconds)
    load(dataImgsFile); %loads predefined cell array called raw_filenames into memory
    num_runs = length(raw_filenames)/expt.numTpPerRun; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
    num_conds_total = size(onsets,2); %onsets will have been loaded into workspace
    thisSubjFname = [thisMvpaDir '/' subjId '_' expt.roiName '_s8mm_wa.mat']; %having this previously saved avoids time-consuming data extraction and preproc
    %if subj structure hasn't been created and saved, do data proc routine (load pattern, detrend, high-pass, filter, z-score)
    if ~exist(thisSubjFname,'file')
        subj = JR_mvpa_load_and_preprocess_raw_data(subjId, expt.name, expt.roiName, expt.roiFname, raw_filenames, num_runs, expt.numTpPerRun);
        save(thisSubjFname, 'subj');
    else
        load(thisSubjFname)
    end
    % initialize all_onsets matrix as conditions x timepoints
    all_onsets = zeros(num_conds_total,size(subj.patterns{end}.mat,2));
    condensed_regs_all = [];
    condensed_runs = [];
    for condNum = 1: num_conds_total
        for trial = 1: length(onsets{condNum})
            % divide by trSecs and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            time_idx = onsets{condNum}(trial)/expt.trSecs+1;
            all_onsets(condNum,time_idx) = 1;
        end
    end
    % condense regressors: all_onsets matrix to one column per trial, instead of per TR (remove TRs that have 0's for all conditions,
    % i.e. rest timepoints)
    trial_counter = 1;
    for i = 1: size(all_onsets,2)
        if ~isempty(find(all_onsets(:,i))) % if not a rest timepoint
            condensed_regs_all(:,trial_counter) = all_onsets(:,i);
            condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
            trial_counter = trial_counter + 1;
        end
    end
    % temporally condense data: select TRs of interest (to correspond with peak post-stim BOLD
    % response) and detect outliers if specified
    all_trials = sum(all_onsets,1); % vector of all trials
    for trNum = 1:length(expt.trWeights)
        data_by_TR(trNum,:,:) = expt.trWeights(trNum)*subj.patterns{end}.mat(:,find(all_trials)+trNum-1);
    end
    temporally_condensed_data = squeeze(sum(data_by_TR(1:trNum,:,:),1));
    clear data_by_TR;
    if flags.remove_artdetect_outliers == 1
        % Exclude trials determined to be outliers by custom ArtDetect script
        % Guide to outlier file cell arrays... Movement thresholds: .2 .25 .3
        % .35 .4 .4 .5 Global signal thresholds: 2 2.5 3 3.5 4 4.5 5
        load([expt.dir '/outlier_indices/' subjId '_outlier_indices']); %load outlier indices
        m_outliers = movement_outlier_trials{flags.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
        gs_outliers = global_signal_outlier_trials{flags.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
        combined_outliers = union(m_outliers,gs_outliers);
        condensed_regs_all(:,combined_outliers) = 0;
        display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
        display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
        display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
    end
    % on-the-fly outlier detection/removal;
    if flags.remove_outlier_trials ~= 0
        mean_across_voxels = mean(temporally_condensed_data,1);
        z_mean_across_voxels = zscore(mean_across_voxels);
        upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
        lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
        all_outliers = union(upper_outliers,lower_outliers)
        condensed_regs_all(:,all_outliers) = 0;
    end
    
    %% assign conditions to train/test classifier on
    condensedRegsOfInterest = tbmvpaGetRegsOfInterest(condensed_regs_all,condNames,expt.studyCondCols);
    %assign conditions to balance classifier, if specified
    trnCondensedRegsToBal = ''; tstCondensedRegsToBal = '';
    if ~isempty(trnSubCondsToBal)
        trnCondensedRegsToBal = tbmvpaGetRegsOfInterest(condensed_regs_all, trnSubCondsToBal, expt.studyCondCols);
        if ~isempty(tstSubCondsToBal)
            tstCondensedRegsToBal = tbmvpaGetRegsOfInterest(condensed_regs_all, tstSubCondsToBal, expt.studyCondCols);
        end
    end
    %xvalBins = tbmvpaSpecifyXvalBins(inputs)
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
    
    %% For num full iterations: use loaded subj patterns, resulting temporally condensed data, and regressors (condensed to one value per trial)
    subj_original = subj; % backup the subj structure before entering loop
    condensed_regs_all_original = condensed_regs_all; % backup condensed_regs_all before entering k-loop
    results = cell(flags.num_full_iter*flags.num_results_iter*flags.num_iter_with_same_data,1); %preallocate results
    currNumTotalIters = 0;  %initialize counter (gets incremented during each classification run-through)
    for fullIterNum = 1:flags.num_full_iter
        subj = subj_original;
        condensed_regs_all = condensed_regs_all_original;
        %% Set up regressors and selectors. perform feature selection
        % initialize regressors object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensedRegsOfInterest);
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
        % "activate" only those trials of interest (from regs_of_interest)
        % before creating cross-validation indices
        active_trials = find(sum(condensedRegsOfInterest)); % find those trials where the sum (across columns) is not zeros
        actives_selector = zeros(1,size(condensed_regs_all,2)); % create new selector object; begin by intializing vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        subj = init_object(subj,'selector','conditions_of_interest'); %initialize new selector object called 'conditions_of_interest'
        subj = set_mat(subj,'selector','conditions_of_interest',actives_selector); % set this selector object to have the values specified by actives_selector
        subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest', 'ignore_jumbled_runs', true,'ignore_runs_zeros',true); % create cross-validation indices (new selector group), using only the the trials specified by actives_selector
        
        % feature selection
        statmap_arg.use_mvpa_ver = true;
        if flags.anova_p_thresh ~= 1
            subj = JR_feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        % run feature selection ANOVA: specify #of voxels (if desired)
        if flags.anova_expt.nVox_thresh ~=0
            subj = JR_feature_select_top_N_vox(subj,'epi_d_hp_z_condensed','conds','runs_xval','expt.nVox_thresh',flags.anova_expt.nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        %% For num results iter: balance, z-score & run classifier for each iteration with same data
        subj_prebalancing = subj; % backup the subj structure before entering loop
        active_trials_prebalancing = active_trials; % backup condensed_regs_all before entering loop
        for resultsIterNum = 1: flags.num_results_iter
            subj = subj_prebalancing;
            active_trials = active_trials_prebalancing;
            new_active_trials =[];
            new_actives_selector= zeros(1,size(condensed_regs_all,2));
            % balance the number of Class A and Class B trials within run,
            % exclude random subsets of trials from each run
            if flags.equate_number_of_trials_in_cond_1_and_2 == 1
                subj = tbm_create_balanced_xvalid_selectors(subj,'conds','runs_xval', trnCondensedRegsToBal, tstCondensedRegsToBal); % creates a new 'runs_xval_bal' selector group
            end
            % we need to get a new list of the 'active trials'
            for runNum = 1:num_runs
                new_active_trials = horzcat(new_active_trials, find(subj.selectors{end-num_runs+runNum}.mat==2));
            end
            new_actives_selector(new_active_trials)=1;
            subj = init_object(subj,'selector','conditions_of_interest_bal_within_runs'); %initialize selector object
            subj = set_mat(subj,'selector','conditions_of_interest_bal_within_runs',new_actives_selector);
            active_trials = new_active_trials;
            
            % z-score again, if specified
            if flags.perform_second_round_of_zscoring == 1  % z-score the data prior to classification (active trials only)
                subj.patterns{5}.mat(:,active_trials) = zscore(subj.patterns{5}.mat(:,active_trials)')';
                display('Performing second round of z-scoring')
            end
            %optimize penalty, if specified
            if flags.optimize_penalty_param == 1 && currNumTotalIters == 0 % find empirically optimal penalty via nested cross-validation
                %run this function during the first pass through; use
                %resulting value for subsequent classifications
                [subj best_penalties penalty_iteration_results] = optimal_pLR_penalty(subj,'epi_d_hp_z_condensed','conds','runs_final_xval','runs_final',classifier_mask,'conditions_of_interest_final','use_iteration_perf',false,'perform_final_classification',false);
                class_args.penalty = best_penalties; % because 'use_iteration_perf' is set to 'false', best_penalties will be only a single value (averaged across the nested cross-validation iterations)
            end
            
            %% for each iter with same data: PERFORM CROSS-VALIDATION (loop relevant if using a stochastic (i.e. non-deterministic) classification algorithm like backpropagation neural nets)
            for iterWithSameDataNum = 1:flags.num_iter_with_same_data
                currNumTotalIters=currNumTotalIters+1; % increment results iteration counter
                if flags.scramble ==1
                    subj = JR_scramble_regressors(subj,'conds','runs','conditions_of_interest_bal_within_runs','conds_scrambled');
                    [subj results{currNumTotalIters}] = cross_validation(subj,'epi_d_hp_z_condensed','conds_scrambled','runs_xval_bal',classifier_mask,class_args);
                else %run the normal classifier without scrambling
                    [subj results{currNumTotalIters}] = cross_validation(subj,'epi_d_hp_z_condensed','conds','runs_xval_bal',classifier_mask,class_args);
                end
                % do some important RAM clean-up
                for runNum = 1:num_runs
                    results{currNumTotalIters}.iterations(runNum).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
                end
                
                
                %% analyze the results in more detail
                %initialize some variables - which are actually stored in
                %results....
                correct_vector = []; desireds_vector = []; guesses_vector = []; prob_class_A_vector = [];
                
                %%this line is a gross hack!!!%%
                if which_traintest >=3 && which_traintest ~= 7
                    num_testing_sets = 1;  %% only summarize and record data from COUNTERMEASURES runs (1st cross-validation testing set)
                else
                    num_testing_sets = num_runs;
                end
                
                for a = 1:num_testing_sets % concatenate the results across all cross-validation testing sets
                    correct_vector = horzcat(correct_vector,results{currNumTotalIters}.iterations(a).perfmet.corrects); % create a vector of the whether or not the classifier's guess on each trial was correct (1=correct; 2=incorrect)
                    desireds_vector = horzcat(desireds_vector,results{currNumTotalIters}.iterations(a).perfmet.desireds); % create a vector of the true class of each trial (i.e., the desired outcome) (1=Class A; 2=Class B)
                    guesses_vector = horzcat(guesses_vector,results{currNumTotalIters}.iterations(a).perfmet.guesses); % create a vector of the classifier's guesses on each trial (1=Class A; 2=Class B)
                    prob_class_A_vector = horzcat(prob_class_A_vector, results{currNumTotalIters}.iterations(a).acts(1,:)); % create a vector of the classifier's scalar probability estimates that each trial is a trial of Class A
                end
                
                overall_accuracy = mean(correct_vector); %accuracy measure based on performance across all trials (rather than averaging the performance of each run)
                overall_hit_rate = mean(correct_vector(desireds_vector==1)); %probability of labeling examples of Cond1 as Cond1
                overall_fa_rate = 1-mean(correct_vector(desireds_vector==2)); %probability of labeling examples of Cond2 as Cond1
                overall_recall = sum(correct_vector(desireds_vector==1))/(sum(correct_vector(desireds_vector==1))+sum(correct_vector(desireds_vector==2)));
                overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate); %measure of classifier sensitivity
                
                % sort by absolute value of classifier "confidence"
                [abs_sorted_diffs abs_ind] = sort(abs(prob_class_A_vector-0.5),2,'descend');
                abs_correct_sorted = correct_vector(abs_ind);
                abs_desireds_sorted = desireds_vector(abs_ind);
                
                
                % compute accuracy for top N % of trials (sorted by
                % classifier 'confidence')
                num_trials(currNumTotalIters) = length(abs_correct_sorted);
                display(['Number of testing trials per bin: ' num2str(num_trials(currNumTotalIters)/2)])
                
                bin_intervals =1:-.05:.05; % twenty bins, ranging from 1.0 (performance across all trials) to 0.5 (performance for Top 5% of trials, ranked based on classifier confidence)
                for acc_bin = 1:20
                    included_trial_inds = 1:ceil(num_trials(currNumTotalIters)*bin_intervals(acc_bin));
                    class_A_inds = find(abs_desireds_sorted==1);
                    class_B_inds = find(abs_desireds_sorted==2);
                    
                    acc_percentiles(acc_bin)= mean(abs_correct_sorted(included_trial_inds));
                end
                
                % print the top 100%, 75%, 50%, and 25% accuracy to command
                % window
                display(acc_percentiles([1 6 11 16]))
                
                % record number of trials in each class at each
                % classification confidence quartile
                class_counts_by_quartile_rank(1,1) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*1.0)))==1);
                class_counts_by_quartile_rank(1,2) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*1.0)))==2);
                class_counts_by_quartile_rank(2,1) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*.75)))==1);
                class_counts_by_quartile_rank(2,2) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*.75)))==2);
                class_counts_by_quartile_rank(3,1) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*.50)))==1);
                class_counts_by_quartile_rank(3,2) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*.50)))==2);
                class_counts_by_quartile_rank(4,1) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*.25)))==1);
                class_counts_by_quartile_rank(4,2) = count(abs_desireds_sorted((1:ceil(num_trials(currNumTotalIters)*.25)))==2);
                
                % sort by signed classifier "confidence" (for ROI curves)
                [sorted_diffs ind] = sort(prob_class_A_vector,2,'descend');
                correct_sorted = correct_vector(ind);
                desireds_sorted = desireds_vector(ind);
                hit_rate = []; fa_rate = [];
                % create continuous ROC function and visualize
                for i = 1:length(sorted_diffs);
                    hit_rate(i) = length(intersect(find(desireds_sorted == 1),[1:i])) / length(find(desireds_sorted == 1));
                    fa_rate(i) = length(intersect(find(desireds_sorted == 2),[1:i])) / length(find(desireds_sorted == 2));
                end
                
                figure(1);
                clf;
                plot(fa_rate,hit_rate,'.-');
                hold on;
                plot([0 1],[0 1],'r');
                xlabel(sprintf('P(%s|%s)',condNames{1},condNames{2}));
                ylabel(sprintf('P(%s|%s)',condNames{1},condNames{1}));
                title([subjId '_' lbl.condStr '_' lbl.dirStr]);
                
                %compute and display area under the ROC curve (AUC); with
                %Hit/FA rates computed across full data continuum
                auc_overall = auroc(hit_rate',fa_rate')
                
                % create ROC function with 80 discrete bins (this allows
                % the ROC curves to be more easily averaged across
                % individuals)
                roc_bin_intervals = .975:-.025:-1;
                for bin_num = 1:80
                    hits_80(bin_num)=length(intersect(find(desireds_sorted == 1),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 1));
                    fas_80(bin_num)=length(intersect(find(desireds_sorted == 2),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 2));
                end
                auc_80_bins = auroc(hits_80',fas_80');
                data_log.overall_acc(currNumTotalIters)=overall_accuracy;
                data_log.hits(currNumTotalIters)=overall_hit_rate;
                data_log.FAs(currNumTotalIters)=overall_fa_rate;
                data_log.d_prime(currNumTotalIters)=overall_d_prime;
                data_log.acc_percentiles(currNumTotalIters,:) = acc_percentiles;
                data_log.penalty_param(currNumTotalIters) = class_args.penalty;
                data_log.class_counts_by_quartile_rank(:,:,currNumTotalIters) = class_counts_by_quartile_rank;
                data_log.auc_overall(currNumTotalIters) = auc_overall;
                data_log.auc_80_bins(currNumTotalIters) = auc_80_bins;
                data_log.roc_80_bin_hits(currNumTotalIters,:)= hits_80;
                data_log.roc_80_bin_fas(currNumTotalIters,:)= fas_80;
            end
            data_log.confmat(subNum,currNumTotalIters,:,:) = multiple_iterations_confusion(results{currNumTotalIters})  ;
        end
        
    end
    
    
    %save outputs as specified
    if flags.save_data_log_as_mat_file ==1;
        %         save_cmd = ['save ' results_data_logs_mat_dir '/' subjId '_' lbl.condStr '.mat data_log flags'];
        %         eval(save_cmd);
        save([lbl.dataLogMatDir '/' subjId '_' lbl.condStr '.mat'], 'data_log', 'flags', 'results');
    end
    if flags.write_data_log_to_text_file==1
        writeDataLogToTextFile(results_data_logs_txt_dir, subjId, condNames, ...
            expt, num_trials, flags, class_args, data_log, currNumTotalIters)
    end
    if flags.generate_importance_maps == 1;
        generateImportanceMaps(class_args, subjId, lbl.impMapDirStr, expt, subj, results, num_testing_sets, condNames, flags, num_runs);
    end
    %record time running
    time2finish = toc/60;
    display(['Finished ' subjId ' in ' num2str(time2finish) ' minutes']);
    %clear all variables except those indicated by keep
    keep  condNames expt lbl flags subNum subj_array trnSubCondsToBal tstSubCondsToBal which_traintest nickname  class_args
end
end




%% Helper function to output importance maps

function [impmap subj] = generateImportanceMaps(class_args, subjId, importance_maps_dir, expt, subj, results, num_testing_sets, condNames, flags, num_runs)
weights = []
impmap = cell(length(condNames),1);
impmap_avg = impmap;
for x = 1:length(results)
    for runNum = 1:num_runs
        switch class_args.train_funct_name
            case 'train_bp'
                weights{x}.iterations(runNum).scratchpad.net.IW{1}  = results{x}.iterations(runNum).scratchpad.net.IW{1};
            case 'train_pLR'
                weights{x}.iterations(runNum).scratchpad.net.IW{1} = results{x}.iterations(runNum).scratchpad.weights';
            case 'train_svdlr'
                weights{x}.iterations(runNum).scratchpad.net.IW{1} = results{x}.iterations(runNum).scratchpad.W';
            otherwise
                weights{x}.iterations(runNum).scratchpad.net.IW{1} = results{x}.iterations(runNum).scratchpad.w';
        end
    end
end

subj = JR_extract_voxel_importance_values_kclass(subj, results,weights);
load([expt.dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
% To create vol_info.mat, run the following command in the matlab workspace:
% vol_info = spm_vol('name_of_image_file_with_same_dimensions_as_data_being_used_by_the_classifier.img'; save vol_info.mat vol_info
for condNum=1:length(condNames)
    impmap{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map i
    if flags.anova_p_thresh == 1 % NO ANOVA VERSION
        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        for testSetNum = 1:num_testing_sets
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            temp{condNum}(voxel_inds)=subj.patterns{end-num_runs+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
            impmap{condNum} = impmap{condNum}+temp{condNum}; %add values cumulatively across iterations
        end
        impmap{condNum} = impmap{condNum}/num_runs*1000; %compute average and multiply by 1000 for scaling
        vol_info.fname = [importance_maps_dir '/' subjId '_' condNames{condNum} '.img'];
        
    else % ANOVA VERSION
        for testSetNum = 1:num_testing_sets
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            voxel_inds{testSetNum} = find(subj.masks{end-num_runs+testSetNum}.mat); %get mask voxel indices
            temp{condNum}(voxel_inds{testSetNum})=subj.patterns{end-num_runs+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
            impmap{condNum} = impmap{condNum}+temp{condNum}; %add values cumulatively across iterations
        end
        %sum across masks to get composite mask (where value of each voxel = number of runs for which that voxel was included)
        composite_mask = zeros(vol_info.dim);
        for maskNum = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
            composite_mask = composite_mask+subj.masks{maskNum}.mat;
        end
        voxels_to_exclude = find(composite_mask<=5);  % exclude voxels that exist for fewer than 5 of the ANOVA masks
        impmap{condNum}(voxels_to_exclude)=0;
        impmap_avg{condNum} = impmap{condNum}./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
        vol_info.fname = [importance_maps_dir '/' subjId '_' condNames{condNum} '_p' num2str(flags.anova_p_thresh) '.img'];
        impmap{condNum} = impmap_avg{condNum};
    end
    spm_write_vol(vol_info,impmap{condNum});
    
    
end
end



%% Helper function to write out datalogs


function [] = writeDataLogToTextFile(results_data_logs_txt_dir, subjId, condNames, ...
    expt, num_trials, flags, class_args, data_log, currNumTotalIters)

filename= [results_data_logs_txt_dir '/' subjId '_' condNames{1} '_vs_' condNames{2} '.txt'];
fid=fopen(filename, 'wt');
fprintf(fid, '%s\r\n', ['subj_id = ' subjId]);
fprintf(fid, '%s\r\n', ['ROI_name = ' expt.roiName]);
fprintf(fid, '%s\r\n', ['expt.dataImgsToUse =' expt.dataImgsToUse]);
fprintf(fid, '%s\r\n', ['TR_weights = ' num2str(expt.trWeights)]);
fprintf(fid, '%s\r\n', ['classification:' condNames{1} ' vs. ' condNames{2}]);
fprintf(fid, '%s\r\n', ['Number of Trials Per Bin: ' num2str(mean(num_trials)/2)]);
fprintf(fid, '%s\r\n', ['flags.perform_second_round_of_zscoring = ' num2str(flags.perform_second_round_of_zscoring)]);
fprintf(fid, '%s\r\n', ['flags.remove_mvpa_outlier_trials (std dev) = ' num2str(flags.remove_outlier_trials)]);
fprintf(fid, '%s\r\n', ['flags.remove_artdetect_outlier_trials (std dev) = ' num2str(flags.remove_artdetect_outliers)]);
fprintf(fid, '%s\r\n', ['flags.artdetect_motion_thresh = ' num2str(flags.artdetect_motion_thresh)]);
fprintf(fid, '%s\r\n', ['flags.artdetect_global_signal_thresh = ' num2str(flags.artdetect_global_signal_thresh)]);

if isfield(class_args, 'penalty')
    fprintf(fid, '%s\r\n', ['penalty param = ' num2str(class_args.penalty)]);
end
fprintf(fid, '\n\n');

for q=1:currNumTotalIters
    fprintf(fid, '%4.4f\t', q);
    fprintf(fid, '%4.4f\t', data_log.auc_overall(q));
    fprintf(fid, '%4.0f\t', (num_trials(q)/2));
    fprintf(fid, '%4.4f\t', data_log.acc_percentiles(q,:));
    fprintf(fid, '\t');
    fprintf(fid, '%4.4f\t', data_log.roc_80_bin_hits(q,:));
    fprintf(fid, '\t');
    fprintf(fid, '%4.4f\t', data_log.roc_80_bin_fas(q,:));
    fprintf(fid, '\n');
end

fprintf(fid, '%s\t', 'mean');
fprintf(fid, '%4.4f\t', mean(data_log.auc_overall));
fprintf(fid, '%4.0f\t', (mean(num_trials(q))/2));
fprintf(fid, '%4.4f\t', mean(data_log.acc_percentiles,1));
fprintf(fid, '\t');
fprintf(fid, '%4.4f\t', mean(data_log.roc_80_bin_hits,1));
fprintf(fid, '\t');
fprintf(fid, '%4.4f\t', mean(data_log.roc_80_bin_fas,1));
fprintf(fid, '\n');
fclose(fid);
end

