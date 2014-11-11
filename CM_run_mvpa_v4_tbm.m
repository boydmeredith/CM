function [res]= CM_run_mvpa_v4_tbm(subj_array, nickname, varargin)
%function [res] = CM_run_mvpa_v4_tbm(subj_array, nickname, varargin)
%inputs:
%subj_array - the subject numbers desired for classification. default is all subs ([1,3:10,12:26])
%nickname - string desired for saving results, if not specified, default is based on input arguments
%varargin - parameters for classification will be automatically determined by CM_mvpa_params if not otherwise specified. Parameters that can be set include: 'condNames', 'which_traintest', 'trWeights', 'roiName', and others
%outputs:
%res - a structure that keeps track of the parameters used for classification as well as the results of each iteration of classifcation

%set default subj_array
if isempty(subj_array)
    subj_array=[1,3:10,12:26];
end

%% for every subject: preproc, classify & output results as specified
for iSub=subj_array
    tic;
    %load relevant info with CM_mvpa_params
    [expt, classArgs, ~,  par] = CM_mvpa_params(iSub, 'ret');
    expt = acceptUserInput(expt, varargin);
    expt.acrossTrs = ismember(0,(expt.trWeights_train==expt.trWeights_test))
    subjId = sprintf(expt.subjIdFormat,iSub);
    expt.thisMvpaDir = fullfile(expt.dir, subjId, expt.mvpaDirStr);
    expt = makeSaveName(expt, subjId, nickname);
    
    % tell the user what is happening
    if subj_array(1) == iSub
        fprintf('\n\n\nPLEASE CHECK VERIFY THESE CLASSIFICATION PARAMETERS!\nContinuing in 10s...\n\n\n');
        display(expt)
        pause(10)
    end
    fprintf('\n\nBeginning subject CM%03d...', iSub);
    
    %% load or create subj and temporally condense. load onsets to use as regressors and condense so that we have one per trial (instead of 1 per TR)
    % mvpa directory within each subject's folder (will be created below if it doesn't exist)
    onsetsFile = fullfile(expt.thisMvpaDir, expt.onsetsFname);
    expt.scanfiles = vertcat(par.swascanfiles.(par.task));
    nScanFiles = length(expt.scanfiles);
    expt.ons = load(onsetsFile); % load in your SPM-formatted onsets file (specifies onsets of each event time in seconds)
    onsets = expt.ons.onsets;
    
    expt.condCols = makeCondCols(expt);
    
    nRuns = nScanFiles/expt.numTpPerRun; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
    nCondsTotal = size(onsets,2);
    
    nConds = length(expt.condNames);
    nExamples = [];
    for i =1:nConds
        nExamples(i) = length(onsets{expt.condCols(i)});
    end
    
    %if subj structure hasn't been created and saved, do data proc routine (load pattern, detrend, high-pass, filter, z-score)
    if ~exist(expt.subjFname,'file')
        subj = CM_mvpa_load_and_preprocess_raw_data(subjId, expt, nRuns, 1)
    else
        load(expt.subjFname)
    end
    nPatts = size(subj.patterns{end}.mat,2)
    
    assert(nScanFiles == nPatts);
    
    % initialize all_onsets matrix as conditions x timepoints
    condOnsets = zeros(nConds,nPatts);
    condensed_runs = [];
    
    for i = 1:nConds %convert from sec to TRs
        condNum = expt.condCols(i);
        time_idx = onsets{condNum}/expt.trSecs+1;% divide by trSecs and add 1 (first timepoint = 0 sec; first TR = 1)
        condOnsets(i,time_idx) = 1;
    end
    all_trials = sum(condOnsets,1);
    restTp = (all_trials==0);    %condense regressors: condOnsets matrix to one column per trial, instead of per TR

    expt.meta_runs = expt.which_traintest(subj.selectors{1}.mat);
    expt.meta_runs(restTp) = 0;
         
    condensedCondRegs = condOnsets(:,expt.meta_runs~=0); %remove TRs that have 0's for all conditions i.e. rest timepoints
    
    % manipulate runs based on which_traintest to set xvalid bins
    condensed_runs = subj.selectors{1}.mat(expt.meta_runs~=0);
    res.subj{iSub}.penalty(1).nVox(1).weights(1).condensed_runs = condensed_runs;
    condensed_runs = expt.which_traintest(condensed_runs);
    
    % temporally condense data: select TRs of interest (to correspond with peak post-stim BOLD
    % response) and detect outliers if specified
    if expt.acrossTrs
        isOneTrainGroupOneTestGroup = (nanmean(ismember(expt.which_traintest,[0 1 2]))==1);
        assert(isOneTrainGroupOneTestGroup);
        expt.meta_runs_train = find(expt.meta_runs==1);
        expt.meta_runs_test = find(expt.meta_runs==2);
    else
       expt.meta_runs_train = find(all_trials);
       expt.meta_runs_test = [];
    end
    
    [temporally_condensed_data] = temporallyCondenseData(subj, expt);
    
    if expt.remove_artdetect_outliers
        % Exclude trials determined to be outliers by custom ArtDetect script
        % Guide to outlier file cell arrays... Movement thresholds: .2 .25 .3
        % .35 .4 .4 .5 Global signal thresholds: 2 2.5 3 3.5 4 4.5 5
        load([expt.dir '/outlier_indices/' subjId '_outlier_indices']); %load outlier indices
        m_outliers = movement_outlier_trials{expt.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
        gs_outliers = global_signal_outlier_trials{expt.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
        combined_outliers = union(m_outliers,gs_outliers);
        condensedCondRegs(:,combined_outliers) = 0;
        display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
        display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
        display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
    end
    
    if expt.remove_outlier_trials ~= 0
        % on-the-fly outlier detection/removal;
        mean_across_voxels = mean(temporally_condensed_data,1);
        z_mean_across_voxels = zscore(mean_across_voxels);
        upper_outliers = find(z_mean_across_voxels> expt.remove_outlier_trials);
        lower_outliers = find(z_mean_across_voxels< -1 * expt.remove_outlier_trials);
        all_outliers = union(upper_outliers,lower_outliers)
        condensedCondRegs(:,all_outliers) = 0;
    end
    
    % assign conditions to train/test classifier on
    %assign conditions to balance classifier, if specified
    trnCondensedRegsToBal = ''; tstCondensedRegsToBal = '';
    if ~isempty(expt.trnSubCondsToBal)
        trnCondensedRegsToBal = tbmvpaGetRegsOfInterest(condensedCondRegs, expt.trnSubCondsToBal, expt.condCols);
        if ~isempty(expt.tstSubCondsToBal)
            tstCondensedRegsToBal = tbmvpaGetRegsOfInterest(condensedCondRegs, expt.tstSubCondsToBal, expt.condCols);
        end
    end
    %xvalBins = tbmvpaSpecifyXvalBins(inputs)
    
    %modify nRuns now that we've altered to make the training and testing
    %sets
    nRuns = length(unique(nonzeros(condensed_runs)));
    
    %% For num full iterations: use loaded subj patterns, resulting temporally condensed data, and regressors (condensed to one value per trial)
    subj_original = subj; % backup the subj structure before entering loop
    condensedCondRegsOrig = condensedCondRegs; % backup condensed_regs_all before entering k-loop
    results = cell(expt.num_full_iter*expt.num_results_iter*expt.num_iter_with_same_data,1); %preallocate results
    currNumTotalIters = 0;  %initialize counter (gets incremented during each classification run-through)
    for fullIterNum = 1:expt.num_full_iter
        %restore from backups
        subj = subj_original;
        condensedCondRegs = condensedCondRegsOrig;
        
        %% Set up regressors and selectors. perform feature selection
        % initialize regressors object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensedCondRegs);
        subj = set_objfield(subj,'regressors','conds','condnames',expt.condNames);
        % add new condensed activation pattern
        subj = duplicate_object(subj,'pattern','epi_d_hp_z','epi_d_hp_z_condensed');
        subj = set_mat(subj,'pattern','epi_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);
        zhist = sprintf('Pattern ''%s'' created by JR custom code','epi_d_hp_z_condensed');
        subj = add_history(subj,'pattern','epi_d_hp_z_condensed',zhist,true);
        % clean up workspace to save RAM
        subj = remove_mat(subj,'pattern','epi_d_hp_z'); % remove original uncondensed pattern (full timeseries)
        
        % update run vector to condensed format.
        subj.selectors{1}.mat = condensed_runs;
        subj.selectors{1}.matsize = size(condensed_runs);
        
        % "activate" only those trials of interest (from regs_of_interest)
        % before creating cross-validation indices
        active_trials = find(sum(condensedCondRegs)); % find those trials where the sum (across columns) is not zeros
        actives_selector = zeros(1,size(condensedCondRegs,2)); % create new selector object; begin by intializing vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        subj = init_object(subj,'selector','conditions_of_interest'); %initialize new selector object called 'conditions_of_interest'
        subj = set_mat(subj,'selector','conditions_of_interest',actives_selector); % set this selector object to have the values specified by actives_selector
        subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest', 'ignore_jumbled_runs', true,'ignore_runs_zeros',true); % create cross-validation indices (new selector group), using only the the trials specified by actives_selector
        
        % feature selection
        statmap_arg.use_mvpa_ver = true;
        if expt.anova_p_thresh ~= 1
            subj = JR_feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',expt.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        if expt.anova_nVox_thresh ~=0
            % run feature selection ANOVA: specify #of voxels (if desired)
            subj = JR_feature_select_top_N_vox(subj,'epi_d_hp_z_condensed','conds','runs_xval','expt.nVox_thresh',expt.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        %% For num results iter: balance, z-score & run classifier for each iteration with same data
        subj_prebalancing = subj; % backup the subj structure before entering loop
        active_trials_prebalancing = active_trials; % backup condensed_regs_all before entering loop
        
        for resultsIterNum = 1: expt.num_results_iter
            
            subj = subj_prebalancing;
            active_trials = active_trials_prebalancing;
            new_active_trials =[];
            new_actives_selector= zeros(1,size(condensedCondRegs,2));
            
            if expt.equate_number_of_trials_in_cond_1_and_2 == 1
                % balance the number of Class 1..N trials within run,
                % exclude random subsets of trials from each run
                subj = tbm_create_balanced_xvalid_selectors(subj,'conds','runs_xval', trnCondensedRegsToBal, tstCondensedRegsToBal); % creates a new 'runs_xval_bal' selector group
            end
            
            % we need to get a new list of the 'active trials'
            for runNum = 1:nRuns
                new_active_trials = horzcat(new_active_trials, find(subj.selectors{end-nRuns+runNum}.mat==2));
            end
            new_actives_selector(new_active_trials)=1;
            subj = init_object(subj,'selector','conditions_of_interest_bal_within_runs'); %initialize selector object
            subj = set_mat(subj,'selector','conditions_of_interest_bal_within_runs',new_actives_selector);
            active_trials = new_active_trials;
            
            if expt.perform_second_round_of_zscoring
                % z-score the data prior to classification (active trials only)
                subj.patterns{5}.mat(:,active_trials) = zscore(subj.patterns{5}.mat(:,active_trials)')';
                display('Performing second round of z-scoring')
            end
            
            if expt.optimize_penalty_param && currNumTotalIters == 0 % find empirically optimal penalty via nested cross-validation
                %run this function during the first pass through; use
                %resulting value for subsequent classifications
                [subj best_penalties penalty_iteration_results] = optimal_pLR_penalty(subj,'epi_d_hp_z_condensed','conds','runs_final_xval','runs_final',classifier_mask,'conditions_of_interest_final','use_iteration_perf',false,'perform_final_classification',false);
                classArgs.penalty = best_penalties; % because 'use_iteration_perf' is set to 'false', best_penalties will be only a single value (averaged across the nested cross-validation iterations)
            end
            
            %% for each iter with same data: PERFORM CROSS-VALIDATION (loop relevant if using a stochastic (i.e. non-deterministic) classification algorithm like backpropagation neural nets)
            
            for iterWithSameDataNum = 1:expt.num_iter_with_same_data
                currNumTotalIters=currNumTotalIters+1; % increment results iteration counter
                
                if expt.scramble
                    subj = JR_scramble_regressors(subj,'conds','runs','conditions_of_interest_bal_within_runs','conds_scrambled');
                    [subj results] = cross_validation(subj,'epi_d_hp_z_condensed','conds_scrambled','runs_xval_bal',classifier_mask,classArgs,'perfmet_functs', expt.perfmetFuncts);
                else %run the normal classifier without scrambling
                    [subj results] = cross_validation(subj,'epi_d_hp_z_condensed','conds','runs_xval_bal',classifier_mask,classArgs,'perfmet_functs', expt.perfmetFuncts);
                end
                
                if expt.generate_importance_maps
                    %store relevant data for creating importance maps later
                    results_IW{currNumTotalIters} = record_weights(results, nRuns, classArgs);
                end
                
                %save results
                res.subj{iSub}.penalty(1).nVox(1).weights(1).iter{currNumTotalIters} = results;
                res.subj{iSub}.penalty(1).nVox(1).weights(1).expt{currNumTotalIters} = expt;
                res.subjArray = subj_array;
                
                
            end
        end
    end
    
    if isempty(expt.xval_iters_for_imp_map)
        expt.xval_iters_for_imp_map = 1:nRuns;
    end
    
    if expt.generate_importance_maps
        if ~exist(expt.impMapDirStr)
            mkdir(expt.impMapDirStr);
        end
        generateImportanceMaps(classArgs, subjId, expt, subj, results,results_IW, nRuns);
    end
    
    %record time running
    time2finish = toc/60;
    display(['Finished ' subjId ' in ' num2str(time2finish) ' minutes']);
    
    clear subj
end
if expt.generate_importance_maps
    write_mean_imp_map(expt);
end
if ~(exist(expt.group_mvpa_dir))
    mkdir(expt.group_mvpa_dir);
end
save(fullfile(expt.group_mvpa_dir, [expt.saveName '.mat']), 'res');

end

function write_mean_imp_map(expt)
cd(expt.impMapDirStr);
for c = 1:length(expt.condNames)
    dfn = dir(fullfile(expt.impMapDirStr, ['*' expt.condNames{c} '*']));
    fn = {dfn.name};
    if length(fn) > length(expt.subj_array), break, end;
    out_name = regexprep(fn{1},'CM\d*','mean');
    spm_imcalc_ui(fn,out_name,'mean(X)',{1,[],[],[]});
    
end
end

function results_IW = record_weights(results, nRuns, classArgs)
for rif = 1:nRuns
    switch classArgs.train_funct_name
        case 'train_pLR'
            this_scratch = results.iterations(rif).scratchpad.weights';
        otherwise
            display('Whoops, don''t know how to archive those feature weights');
    end
    results_IW.iterations(rif).scratchpad.net.IW{1} = this_scratch;
end
end

function expt = acceptUserInput(expt,  args)
%% allows user to make changes from the parameters specified in CM_mvpa_params without having to change everything every time
p = inputParser;
p.addParamValue('roiName', expt.roiName, @(x) ischar(x));
p.addParamValue('trWeights_train', expt.trWeights_train, @(x) isnumeric(x));
p.addParamValue('trWeights_test', expt.trWeights_test, @(x) isnumeric(x));
p.addParamValue('which_traintest', expt.which_traintest, @(x) isnumeric(x));
p.addParamValue('condNames', expt.condNames, @(x) iscell(x));
p.addParamValue('trnSubCondsToBal', expt.trnSubCondsToBal, @(x) iscell(x) || ischar(x));
p.addParamValue('tstSubCondsToBal', expt.tstSubCondsToBal, @(x) iscell(x) || ischar(x));
p.addParamValue('num_results_iter',expt.num_results_iter, @(x) isnumeric(x));
p.addParamValue('generate_importance_maps',expt.generate_importance_maps, @(x) (x==1 || x ==0));
p.addParamValue('xval_iters_for_imp_map',[],@(x) isnumeric(x));
p.parse(args{:});
res = p.Results;
expt = mergestructs(res, expt);
end


%% sets the indices of the conditions that need to be examined
function condCols = makeCondCols(expt)
condCols = zeros(1,length(expt.condNames));
for i = 1:length(condCols)
    thisName = expt.condNames{i};
    condCols(i) = find(ismember(expt.ons.names,thisName));
end
end

%% Helper function to output importance maps

function [impmap subj] = generateImportanceMaps(class_args, subjId,  expt, subj, results, results_IW,  numRuns)

assert(mean(ismember(expt.xval_iters_for_imp_map,1:numRuns))==1); %check xval_iters_for_imp_map are all real runs
impmap = cell(length(expt.condNames),1);
impmap_avg = impmap;

load([expt.dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
subj = JR_extract_voxel_importance_values_kclass(subj, results,results_IW);

% To create vol_info.mat, run the following command in the matlab workspace:
% vol_info = spm_vol('name_of_image_file_with_same_dimensions_as_data_being_used_by_the_classifier.img'; save vol_info.mat vol_info
for condNum=1:length(expt.condNames)
    impmap{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map i
    if expt.anova_p_thresh == 1 % NO ANOVA VERSION
        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        for testSetNum = 1:expt.xval_iters_for_imp_map
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            temp{condNum}(voxel_inds)=subj.patterns{end-numRuns+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
            impmap{condNum} = impmap{condNum}+temp{condNum}; %add values cumulatively across iterations
        end
        impmap{condNum} = impmap{condNum}/numRuns*1000; %compute average and multiply by 1000 for scaling
        vol_info.fname = [expt.impMapDirStr '/' subjId '_' expt.condNames{condNum} '_xals' strrep(num2str(expt.xval_iters_for_imp_map),' ','') '.img'];
        
    else % ANOVA VERSION
        for testSetNum = 1:expt.xval_iters_for_imp_map
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            voxel_inds{testSetNum} = find(subj.masks{end-numRuns+testSetNum}.mat); %get mask voxel indices
            temp{condNum}(voxel_inds{testSetNum})=subj.patterns{end-numRuns+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
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
        vol_info.fname = [expt.impMapDirStr '/' subjId '_' expt.condNames{condNum} '_p' num2str(expt.anova_p_thresh) '.img'];
        impmap{condNum} = impmap_avg{condNum};
    end
    spm_write_vol(vol_info,impmap{condNum});
    
    
end





end


function expt = makeSaveName(expt, subjId, nickname)

if isempty(nickname)
    condNamesStr = strrep(strjoin(expt.condNames),',','V');
    trTeStr = strrep(num2str(expt.which_traintest),'  ','_');
    trWStr_tr = strrep(num2str(expt.trWeights_train),' ','_');
    trWStr_te = strrep(num2str(expt.trWeights_test),' ','_');
    trWStr = ['tr' trWStr_tr '_te' trWStr_te];
    roiStr = regexprep(expt.roiName,'(\w*).*','$1');
    nickname = sprintf('conds_%s_trTe%s_trW%s_roi%s',condNamesStr,trTeStr,trWStr,roiStr);
    
end
expt.saveName = [expt.saveName nickname];
expt.impMapDirStr=[expt.dir '/mvpa_results/importance_maps/' expt.saveName];
expt.roiFname = [expt.dir '/Masks/' expt.roiName]; % specific path to mask file
expt.subjFname = [expt.thisMvpaDir '/' subjId '_' expt.roiName '_s8mm_wa.mat']; %having this previously saved avoids time-consuming data extraction and preproc

end


function temporally_condensed_data = temporallyCondenseData(subj, expt)
data_by_TR_train=[];
data_by_TR_test=[];

for iTr = expt.trsToAverageOver
    data_by_TR_train(iTr,:,:) = expt.trWeights_train(iTr)*subj.patterns{end}.mat(:,expt.meta_runs_train+iTr-1);
    data_by_TR_test(iTr,:,:) = expt.trWeights_test(iTr)*subj.patterns{end}.mat(:,expt.meta_runs_test+iTr-1);
    
end

temporally_condensed_data_train = squeeze(sum(data_by_TR_train(expt.trsToAverageOver,:,:),1));
clear data_by_TR_train;
temporally_condensed_data_test = squeeze(sum(data_by_TR_test(expt.trsToAverageOver,:,:),1));
clear data_by_TR_test;
temporally_condensed_data = horzcat(temporally_condensed_data_train, temporally_condensed_data_test);

end

