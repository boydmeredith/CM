function [results]= CM_run_mvpa_across_subjects(subj_array, trte, varargin)
%function [res] = CM_run_mvpa_across_subjects(subj_array, trte, varargin)
%inputs: subj_array - the subject numbers desired for classification. default
%is all subs ([1,3:10,12:26]) varargin - parameters for classification will be
%automatically determined by CM_mvpa_params if not otherwise specified.
%Parameters that can be set include: 'condNames', 'which_traintest',
%'trWeights', 'roiName', and others outputs: res - a structure that keeps track
%of the parameters used for classification as well as the results of each
%iteration of classifcation

%determine which classifier to use
class_args.train_funct_name = 'train_pLR'; 
class_args.test_funct_name = 'test_pLR'; 
class_args.penalty = 10000;
    
%set default subj_array
if isempty(subj_array) 
    subj_array=[1,3:10,12:26]; 
end
expt.subj_array =subj_array;

%initialize concatenated subjects variable
concatAllSubs.pattern = []; concatAllSubs.regressors = [];
concatAllSubs.actives = []; concatAllSubs.subjIdx = [];
concatAllSubs.condensed_runs = [];

% for each subject, set up struct and put into voxel ind
for subIdx=1:length(subj_array) tic; subnum = subj_array(subIdx);
    %get expt structure, to determine classification scheme
    [expt classargs s par] = CM_mvpa_params(subnum, 'ret');
    %incorporate user specified variables into default expt structure
    expt = accept_user_input(expt, varargin);
    %handle naming relevant files
    nickname = make_savename(expt);
    expt.saveName = [expt.saveName nickname];
    expt.impmapdirstr=[expt.dir '/mvpa_results/within_subj_class/importance_maps/' expt.saveName];
    expt.roiFname = [expt.dir '/masks/' expt.roiName];
    subjid = sprintf(expt.subjIdFormat,subnum);
    thismvpadir = fullfile(expt.dir, subjid, expt.mvpaDirStr);
    onsetsfile = fullfile(thismvpadir, expt.onsetsFname);
    expt.scanfiles = vertcat(par.swascanfiles.(par.task));
    nscanfiles = length(expt.scanfiles);
    nruns = nscanfiles/expt.numTpPerRun;     
    % dataimgsfile = fullfile(thismvpadir, expt.dataimgstouse);
    %load(dataimgsfile); %loads predefined cell array called expt.scanfiles into memory
    load(onsetsfile); % load in your spm-formatted onsets file (specifies onsets of each event time in seconds)
    
    %if subj structure hasn't been created and saved, do data proc routine (load pattern, detrend, high-pass, filter, z-score)
    expt.subjFname = [thismvpadir '/' subjid '_' expt.roiName '_s8mm_wa.mat']; %having this previously saved avoids time-consuming data extraction and preproc
    if ~exist(expt.subjFname,'file'), 
        subj = CM_mvpa_load_and_preprocess_raw_data(subjid, expt, nruns, 1)
    else
        load(expt.subjFname)
    end
       
    expt.condCols = make_cond_cols(expt, names);
    % extract info about conditions from onsets file (uses spm-based coding
    % of the onsets (in seconds) of trials in each condition of interest)
    num_conds = size(onsets,2);
    all_regs = zeros(num_conds,nruns*expt.numTpPerRun); % initialize regs matrix as conditions x timepoints
    
    for cond = 1: num_conds %(exclude last condition ("no_response" trials)
        for trial = 1: length(onsets{cond})
            time_idx = round(onsets{cond}(trial)/2)+1; % divide by 2 and add 1 to convert back from sec to trs (first timepoint = 0 sec; first tr = 1)
            all_regs(cond,time_idx) = 1;
        end
    end
    
    % condense all_onsets matrix from having one column per TR to having
    % once column per trial (remove TRs that have 0's for all
    % conditions, i.e. rest timepoints)
    condensed_regs_all = [];
    condensed_runs = [];
    trial_counter = 1;
    for i = 1: size(all_regs,2)
        if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
            %             condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
            condensed_regs_all(:,trial_counter) = all_regs(:,i);
            condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
            trial_counter = trial_counter + 1;
        end
    end
    condensed_runs = expt.which_traintest(condensed_runs);
    
    condensed_regs_of_interest = [];
    for i=1:length(expt.condCols)
        condensed_regs_of_interest(i,:) = condensed_regs_all(expt.condCols(i),:);
    end  
    
    % temporally condense data: select TRs of interest (to correspond with peak post-stim BOLD
    % response) and detect outliers if specified
    all_trials = sum(all_regs,1); % vector of all trials
    for trNum = expt.trsToAverageOver
        data_by_TR(trNum,:,:) = expt.trWeights(trNum)*subj.patterns{end}.mat(:,find(all_trials)+trNum-1);
    end
    temporally_condensed_data = squeeze(sum(data_by_TR(expt.trsToAverageOver,:,:),1));
    clear data_by_TR;
    
    %% look for outliers
    if expt.remove_artdetect_outliers
        % Exclude trials determined to be outliers by custom ArtDetect script
        % Guide to outlier file cell arrays... Movement thresholds: .2 .25 .3
        % .35 .4 .4 .5 Global signal thresholds: 2 2.5 3 3.5 4 4.5 5
        load([expt.dir '/outlier_indices/' subjId '_outlier_indices']); %load outlier indices
        m_outliers = movement_outlier_trials{expt.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
        gs_outliers = global_signal_outlier_trials{expt.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
        combined_outliers = union(m_outliers,gs_outliers);
        condensed_regs_of_interest(:,combined_outliers) = 0;
        display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
        display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
        display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
    end
    
    
    
    if expt.equate_number_of_trials_in_cond_1_and_2 == 1
        % assign conditions to train/test classifier on
        %assign conditions to balance classifier, if specified
        for i = unique(condensed_runs)
            ix = find(condensed_runs==i);
            condcounts = sum(condensed_regs_of_interest(:,ix),2);
            min_counts = min(condcounts)
            for c = 1:size(condensed_regs_of_interest,1)
                % indices of active timepoints for this condition
                cur_cond_idx = find(condensed_regs_of_interest(c,:) & condensed_runs==i);
                % how many surplus timepoints in this condition
                cur_surplus = length(cur_cond_idx) - min_counts;
                % SAMPLE can't return 0 samples, so skip to the next condition if
                % there are no surplus timepoints
                if ~cur_surplus, continue, end
                % pick a random subset to remove, bringing the number of active
                % timepoints for this condition down to MIN_COUNTS
                %remove_idx = sample(cur_cond_idx,cur_surplus);  % JR commmented out
                random_ints = randperm(length(cur_cond_idx)); %JR added
                remove_idx = cur_cond_idx(random_ints(1:cur_surplus));%JR added
                condensed_regs_of_interest(c,remove_idx) = 0;
            end % c nConds
            
            for c=1:size(condensed_regs_of_interest,1)
                % now, all the conditions should have the same number of timepoints
                assert(count(condensed_regs_of_interest(c,ix))==min_counts);
            end % c nConds
            
            
        end
        % re-confirm that none of the timepoints are overactive
        [isbool isrest isover] = check_1ofn_regressors(condensed_regs_of_interest);
        assert(~isover)
        
    end
      
    %% Set up regressors and selectors. perform feature selection
    % initialize regressors object
    subj = init_object(subj,'regressors','conds');
    subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
    subj = set_objfield(subj,'regressors','conds','condnames',expt.condNames);
    
    % add new condensed activation pattern
    subj = duplicate_object(subj,'pattern','epi_d_hp_z','epi_d_hp_z_condensed');
    subj = set_mat(subj,'pattern','epi_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);
    zhist = sprintf('Pattern ''%s'' created by JR custom code','epi_d_hp_z_condensed');
    subj = add_history(subj,'pattern','epi_d_hp_z_condensed',zhist,true);
    
    % clean up workspace to save RAM
    subj = remove_mat(subj,'pattern','epi_d_hp_z'); % remove original uncondensed pattern (full timeseries)
    clear mean_data;
    
    % update run vector to condensed format
    subj.selectors{1}.mat = condensed_runs;
    subj.selectors{1}.matsize = size(condensed_runs);
    
    % "activate" only those trials of interest (from regs_of_interest)
    % before creating cross-validation indices
    nonactive_trials = find(sum(condensed_regs_of_interest)==0);
    active_trials = find(sum(condensed_regs_of_interest)); % find those trials where the sum (across columns) is not zeros
    actives_selector = zeros(1,size(condensed_regs_of_interest,2)); % create new selector object; begin by intializing vector of all zeros
    actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
    actives_selector(nonactive_trials) = 0; % remove all non-"regs_of_interst" trials (set to zero)
    subj = init_object(subj,'selector','conditions_of_interest'); %initialize new selector object called 'conditions_of_interest'
    subj = set_mat(subj,'selector','conditions_of_interest',actives_selector); % set this selector object to have the values specified by actives_selector
    %subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest', 'ignore_jumbled_runs', true,'ignore_runs_zeros',true); % create cross-validation indices (new selector group), using only the the trials specified by actives_selector
       
    if expt.perform_second_round_of_zscoring
        % z-score the data prior to classification (active trials only)
        subj.patterns{5}.mat(:,active_trials) = zscore(subj.patterns{5}.mat(:,active_trials)')';
        display('Performing second round of z-scoring')
    end    

    %% put into subj-general subj struct
    voxel_inds{subIdx} = find(subj.masks{end}.mat);
    if subIdx == 1
        common_voxel_inds = voxel_inds{1};
        voxels_to_include_current = find(ismember(common_voxel_inds,voxel_inds{subIdx}));
        voxels_to_include_archive = find(ismember(common_voxel_inds,voxel_inds{subIdx}));
        concatAllSubs.pattern = subj.patterns{end}.mat(voxels_to_include_current,active_trials);
    else
        % keep record common voxel inds up to this point
        % find which voxels of the continuously incrementing array exists for the current subject
        % find which of the present subject's voxels are in common with the previous collection
        % find which of the previous collection of voxels are in common with the current conjunction
        archive_voxel_inds = common_voxel_inds;
        common_voxel_inds = intersect(common_voxel_inds,voxel_inds{subIdx});
        voxels_to_include_current = find(ismember(voxel_inds{subIdx},common_voxel_inds));
        voxels_to_include_archive = find(ismember(archive_voxel_inds,common_voxel_inds));
        
        display([num2str(length(common_voxel_inds)) ' common voxels remaining in the mask'])
        concatAllSubs.pattern = horzcat(concatAllSubs.pattern(voxels_to_include_archive,:),subj.patterns{end}.mat(voxels_to_include_current,active_trials));
    end
    
    
    concatAllSubs.subjIdx = horzcat(concatAllSubs.subjIdx, repmat(subIdx,1,length(active_trials)));
    concatAllSubs.regressors = horzcat(concatAllSubs.regressors, subj.regressors{1}.mat(:,find(subj.selectors{2}.mat)));  %get regressors only for active trials
    %concatAllSubs.actives = horzcat(concatAllSubs.actives, subj.selectors{2}.mat);
    concatAllSubs.condensed_runs = horzcat(concatAllSubs.condensed_runs, condensed_runs(:,find(subj.selectors{2}.mat)));

    number_of_trials(subIdx)= length(subj.selectors{2}.mat);
    
    clear condensed_regs_of_interest condensed_regs_all condensed_runs
    
end

subj.patterns{end}.mat = concatAllSubs.pattern;
subj.patterns{end}.matsize = size(concatAllSubs.pattern);
clear concatAllSubs.pattern

subj.regressors{1}.mat = concatAllSubs.regressors;
subj.regressors{1}.matsize = size(concatAllSubs.regressors);
clear concatAllSubs.regressors

if strcmp(trte,'EXCM')
    expt.saveName = ['trEXteCM_' expt.saveName];
    subj.selectors{1}.mat = ismember(concatAllSubs.condensed_runs, 5:8)+1;
    subj.selectors{1}.matsize = size(concatAllSubs.condensed_runs);
elseif strfind(trte,'minus')
    str_ix = strfind(trte,'minus')+length('minus');	
    s_to_exclude = str2num(trte(str_ix:str_ix+1));
    s_to_exclude_ix = find(subj_array==s_to_exclude); %convert to index
    subj.selectors{1}.mat = ones(size(concatAllSubs.subjIdx));
    subj.selectors{1}.mat(concatAllSubs.subjIdx==s_to_exclude_ix) = 2;		
    if strncmp(trte,'minus',5)
    	expt.saveName = ['trAllteAll_' expt.saveName];
    elseif strncmp(trte, 'EXminus',7)
    	expt.saveName = ['trEXteAll_' expt.saveName];
    	subj.selectors{1}.mat(ismember(concatAllSubs.condensed_runs, 5:8) & concatAllSubs.subjIdx ~= s_to_exclude_ix)=0;
    elseif strncmp(trte, 'CMminus',7)
    	expt.saveName = ['trCMteAll_' expt.saveName];
    	subj.selectors{1}.mat(ismember(concatAllSubs.condensed_runs, 1:4) & concatAllSubs.subjIdx ~= s_to_exclude_ix)=0;
    end
else
    subj.selectors{1}.mat = concatAllSubs.subjIdx;
end
clf;
imagesc([concatAllSubs.subjIdx; concatAllSubs.regressors; concatAllSubs.condensed_runs;subj.selectors{1}.mat])
hold on; colorbar; title('runs selector');

%subj = remove_mat(subj,'selectors','conditions_of_interest'); % remove original uncondensed pattern (full timeseries)
subj = create_xvalid_indices(subj,'runs', 'ignore_runs_zeros', true);

subj.masks{1}.mat = zeros(subj.masks{1}.matsize);
subj.masks{1}.mat(common_voxel_inds) = 1;

%% Classification
% run feature selection ANOVA: specify pvalue (if desired)
statmap_arg.use_mvpa_ver = true;
if expt.anova_p_thresh ~= 1
    
    subj = JR_feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',expt.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
    %subj = JR_feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',expt.anova_p_thresh, 'greater_than', true,'statmap_funct','statmap_RFE');
    classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
else
    classifier_mask = subj.masks{1}.name; % use original mask
end

% run feature selection ANOVA: specify #of voxels (if desired)
if expt.anova_nVox_thresh ~=0
    
    subj = JR_feature_select_top_N_vox(subj,'epi_d_hp_z_condensed','conds','runs_xval','nVox_thresh',expt.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
    %subj = JR_feature_select_top_N_vox_iterative(subj,'epi_d_hp_z_condensed','conds','runs_xval','nVox_thresh',expt.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
    classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
else
    classifier_mask = subj.masks{1}.name; % use original mask
end

subj.patterns{end}.mat = single(subj.patterns{end}.mat);  % make pattern a 'single' matrix to save ram and speed up classification (Note: doesn't work with backprop, though)

%run cross-validated classification
[subj results{1}] = cross_validation(subj,'epi_d_hp_z_condensed','conds','runs_xval',classifier_mask,class_args);

%save condensed runs for identification of ex and cm trials
results{1}.condensed_runs = concatAllSubs.condensed_runs;

%save results
save (fullfile(expt.group_mvpa_dir, expt.saveName), 'results');

end



%% allows user to make changes from the parameters specified in CM_mvpa_params without having to change everything every time
function expt = accept_user_input(expt,  args)
p = inputParser;
p.addParamValue('roiName', expt.roiName, @(x) ischar(x));
p.addParamValue('trWeights', expt.trWeights, @(x) isnumeric(x));
p.addParamValue('which_traintest', expt.which_traintest, @(x) isnumeric(x));
p.addParamValue('condNames', expt.condNames, @(x) iscell(x));
p.addParamValue('trnSubCondsToBal', expt.trnSubCondsToBal, @(x) iscell(x) || ischar(x));
p.addParamValue('tstSubCondsToBal', expt.tstSubCondsToBal, @(x) iscell(x) || ischar(x));
p.addParamValue('num_results_iter',expt.num_results_iter, @(x) isnumeric(x));
p.parse(args{:});
res = p.Results;
expt = mergestructs(res, expt);
end


%% sets the indices of the conditions that need to be examined
function condCols = make_cond_cols(expt, names)
condCols = zeros(1,length(expt.condNames));
for i = 1:length(condCols)
    thisName = expt.condNames{i};
    condCols(i) = find(ismember(names,thisName));
end
end

%% Helper function to output importance maps

function [impmap subj] = generateImportanceMaps(class_args, subjId, importance_maps_dir, expt, subj, results, num_testing_sets, numRuns)
weights = []
impmap = cell(length(expt.condNames),1);
impmap_avg = impmap;
for x = 1:length(results)
    for runNum = 1:numRuns
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
for condNum=1:length(expt.condNames)
    impmap{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map i
    if expt.anova_p_thresh == 1 % NO ANOVA VERSION
        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        for testSetNum = 1:num_testing_sets
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            temp{condNum}(voxel_inds)=subj.patterns{end-numRuns+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
            impmap{condNum} = impmap{condNum}+temp{condNum}; %add values cumulatively across iterations
        end
        impmap{condNum} = impmap{condNum}/numRuns*1000; %compute average and multiply by 1000 for scaling
        vol_info.fname = [importance_maps_dir '/' subjId '_' expt.condNames{condNum} '.img'];
        
    else % ANOVA VERSION
        for testSetNum = 1:num_testing_sets
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
        vol_info.fname = [importance_maps_dir '/' subjId '_' expt.condNames{condNum} '_p' num2str(expt.anova_p_thresh) '.img'];
        impmap{condNum} = impmap_avg{condNum};
    end
    spm_write_vol(vol_info,impmap{condNum});
    
    
end
end

function savename = make_savename(expt)
condNamesStr = strrep(strjoin(expt.condNames),',','V');
trTeStr = strrep(num2str(expt.which_traintest),'  ','_');
trWStr = strrep(num2str(expt.trWeights),' ','_');
savename = sprintf('acrossSubs_conds_%s_trTe%s_trW%s_roi%s',condNamesStr,trTeStr,trWStr,expt.roiName);


end

