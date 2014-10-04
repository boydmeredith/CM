function [expt, classArgs, S, par]= CM_mvpa_params(subj_id, task)

% establish parameters for mvpa analysis
% <subj_id> - identifier for the given subject. Can be numerical or a
% string
% <task> 'perc' or 'mnem'

%% EXPT SPECIFIC INFO
expt.name = 'CM_ret';
expt.dir = '/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional/'; % top-level expt. directory where individual subject folders live
expt.mvpaDirStr = 'mvpa';
%expt.dataImgsToUse = 'raw_filenames.mat'; % .mat file containing names of all functional images (must exist for each subject; can be created by running cellstr(SPM.xY.P) on subject's SPM.mat file)
expt.numTpPerRun = 256; % number of TRs per scanning run (coded for fixed value; adjust code structure if variable TR counts across runs)
%expt.roiName = 'rLR_PrePostCentralandSMA';  % name of mask to be used for voxel selection (can be a small ROI, a whole-brain mask, or anywhere in between)
%flags determine how classifier wrapper behaves
expt.num_full_iter = 1; % number of times to run the entire classification process, including feature selection
expt.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
expt.num_results_iter = 10; % number of times to run the post-feature selection classification process (select subset of the data and train/test classifier)
expt.equate_number_of_trials_in_cond_1_and_2 = 1; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
expt.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
expt.anova_nVox_thresh = 0; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
expt.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification (because trial selection and TR selection undoes initial z-scoring)
expt.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
expt.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
expt.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
expt.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
expt.generate_importance_maps = 0; % 1=generate importance maps based on classification weights (scaled by the mean univariate activity of each condition)
expt.write_data_log_to_text_file=0; % save a .txt file that logs critical classification performance data
expt.save_data_log_as_mat_file =1; % save a .mat file that logs critical classification performance data
expt.optimize_penalty_param = 0; % 1 = peformed nested cross-validation for empirical penalty optimization
expt.scramble = 0;
%class args determine which classifier to use and what it should do
plrArgs.train_funct_name = 'train_pLR';
plrArgs.test_funct_name = 'test_pLR';
plrArgs.penalty = 10;
linsvmArgs.train_funct_name = 'train_svm';
linsvmArgs.test_funct_name = 'test_svm';
% linsvmArgs.kernel_type = 0; %0=linear, 1=polynomial, 2=rbf, 4=sigmoid
polysvmArgs.train_funct_name = 'train_svm';
polysvmArgs.test_funct_name = 'test_svm';
polysvmArgs.kernel_type = 1;
rbfsvmArgs.train_funct_name = 'train_svm';
rbfsvmArgs.test_funct_name = 'test_svm';
rbfsvmArgs.kernel_type =2;
bpArgs.train_funct_name = 'train_bp';
bpArgs.test_funct_name = 'test_bp';
bpArgs.nHidden = 1;
AA_subjs = {'CM004', 'CM008', 'CM010', 'CM013', 'CM019', 'CM020', 'CM021', 'CM025'};
all_subjs = {'CM001','CM003','CM004','CM005','CM006','CM007','CM008','CM009','CM010','CM012','CM013','CM014','CM015','CM016','CM017','CM018','CM019','CM020','CM021','CM022','CM023','CM024','CM025','CM026'};
highEXtoEXClassPerfSubs = {'CM001' ,   'CM005' ,   'CM006'  ,  'CM009' ,   'CM010'   , 'CM012'   , 'CM014' , 'CM015' ,   'CM019' ,   'CM022'    ,'CM023' ,   'CM025'};
highClassPerfSubIdx = [1     4     5     8     9    10    12    13    17    20    21    23];
expt.subj_array = [1,3:10,12:26];
classArgs = plrArgs;
% do some string handling for directory organization
% expt.condStr = [repmat(sprintf(['%s', '_vs_'], expt.condNames{1:end-1}),1, ~isscalar(expt.condNames)), sprintf('%s', expt.condNames{end})];
% expt.trStr = ['TRs' strrep(num2str(expt.trWeights),' ','')];
% expt.penStr = '';
% if isfield(classArgs,'penalty')
%     if expt.optimize_penalty_param == 1
%         expt.penStr = 'OPTIMAL_pen_';
%     else
%         expt.penStr = ['_pen' num2str(classArgs.penalty) '_'];
%     end
% end
expt.saveName = '';
if expt.scramble ==1
    expt.saveName = 'SCRAMBLED_'
end

%expt.roiName = 'rLR_Precentral_aal';
expt.roiName = 'SEPT09_MVPA_MASK_resliced4mm';
expt.trnSubCondsToBal = {};
expt.tstSubCondsToBal = {};

expt.trsToAverageOver = [1 2 3 4 5 6];
expt.trWeights = [0 0 .33 .34 .33 0];
expt.nVox_thresh = 0;
expt.subjIdFormat = 'CM%03d';
expt.onsetsFname = 'mvpa_ons.mat';
expt.trSecs = 2; %used to convert from seconds to TRs
expt.perfmetFuncts = 'perfmet_maxclass';
expt.group_mvpa_dir = [expt.dir '/mvpa_results'];
expt.which_traintest = [1 2 3 4 0 0 0 0];
expt.condNames = {'hit', 'cr'};
expt.subjFname = fullfile(expt.dir, subj_id, expt.mvpaDirStr,[subj_id '_' expt.roiName '.mat']); %having this previously saved avoids time-consuming data extraction and preproc

%% establish general parameters
S.idxTr = [];
S.idxTe = [];
par = CM_Params(subj_id, task);
S.exp_name = 'CM';

%% directories
S.subj_id = par.substr;
S.expt_dir = par.funcdir;
S.mvpa_dir = [par.subdir '/mvpa2'];
S.anat_dir = [par.subdir '/anatomy'];
S.importance_maps_dir=[S.expt_dir 'mvpa_results/ImpMaps_' date  ];
S.group_mvpa_dir = [S.expt_dir 'mvpa_results2'];


S.workspace_dir = [par.subdir '/' 'mvpa_workspace'];

%% preprocessing
S.preprocType = 'spm'; % 'spm' for spm preprocessing, 'knk' for kendrick preprocessing
S.loadAndPreprocFunctName = 'JR_mvpa_load_and_preprocess_raw_data';
%% tasks
S.trainTask = 'exLeftHand';
S.testTask = 'exLeftHand';

%% cross-validation scheme
if strcmp(S.trainTask, S.testTask)
    S.xval = 1;
    S.thisSelector =  'randomNFold_xval'; % cross validation. 
else
    S.xval = 0;
    S.thisSelector = 'TrainTestOneIterGroup'; % train on one group, test on another
end

%% information specific to the training and testing of specific sets of data

% parTr = par for the training data
% S.onsetsTrainDir = location of training onsets
% S.condsTrain = conditions on which to train
% S.dnCondsTrain = conditions which which to denoise, if denoising is used
% S.TrainRuns = runs of data on which to train
% S.durTrain = duration of training
% S.filenames_train = names of data images to use for training
% S.idxTr = behavioral information for training task
    
if strcmp(S.trainTask,'exLeftRight')
    parTr = CM_Params(subj_id, 'ex');
    S.onsetsTrainDir = [S.mvpa_dir];
    S.condsTrain = {{'exLeftHand'}  {'exRightHand'}} ;
% %     S.dnCondsTrain = {{'face'}  {'house'} {'noise'}};
    S.TrainRuns = par.scansSelect.(par.task);
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swascanfilesByRun.(par.task){S.TrainRuns});
% %     %S.filenames_train = vertcat(par.betasByRun{S.TrainRuns});
    [~, ~,S.idxTr] = CM_fMRIBehAnalysis(par, 'ex');
elseif strcmp(S.trainTask,'exLeftHand')
    parTr = CM_Params(subj_id, 'ex');
    S.onsetsTrainDir = [S.mvpa_dir];
    S.condsTrain = {{'exLeftIndex'}  {'exLeftMiddle'}} ;
% %     S.dnCondsTrain = {{'face'}  {'house'} {'noise'}};
    S.TrainRuns = par.scansSelect.(par.task);
    S.durTrain = sum(par.(par.task).numvols(S.TrainRuns)) * par.TR;
    S.filenames_train = vertcat(par.swascanfilesByRun.(par.task){S.TrainRuns});
% %     %S.filenames_train = vertcat(par.betasByRun{S.TrainRuns});
    [~, ~,S.idxTr] = CM_fMRIBehAnalysis(par, 'ex');
end
%% testing
if strcmp(S.testTask,'cmLeftRight')
    S.onsetsTestDir = [S.mvpa_dir];
    S.condsTest = {{'cmLeftHand'}  {'cmRightHand'}} ;
    S.TestRuns = par.scansSelect.(par.task);
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.swascanfilesByRun.(par.task){S.TrainRuns});
% %     S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    [~, ~,S.idxTe] = CM_fMRIBehAnalysis(par, 'cm');
elseif strcmp(S.testTask,'cmLeftHand')
    S.onsetsTestDir = [S.mvpa_dir];
    S.condsTest = {{'cmLeftIndex'}  {'cmLeftMiddle'}} ;
    S.TestRuns = par.scansSelect.(par.task);
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.swascanfilesByRun.(par.task){S.TrainRuns});
% %     S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    [~, ~,S.idxTe] = CM_fMRIBehAnalysis(par, 'cm');
elseif strcmp(S.testTask,'exLeftHand')
    S.onsetsTestDir = [S.mvpa_dir];
    S.condsTest = {{'exLeftIndex'}  {'exLeftMiddle'}} ;
    S.TestRuns = par.scansSelect.(par.task);
    S.durTest = sum(par.(par.task).numvols(S.TestRuns)) * par.TR;
    S.filenames_test = vertcat(par.swascanfilesByRun.(par.task){S.TrainRuns});
% %     S.filenames_test = vertcat(par.betasByRun{S.TestRuns});
    [~, ~,S.idxTe] = CM_fMRIBehAnalysis(par, 'ex');

end

S.condnames = S.condsTrain;
S.regName = 'conds';


%% Smoothing Parameters
S.funcType = 2;
S.smoothTxt = { 'unsmoothed' 'smoothed' 'native'};
switch S.funcType
    case 1
    par.filesForPatterns = par.wascanfiles;
    case 2
    par.filesForPatterns = par.swascanfiles;
    case 3
    par.filesForPatterns = par.ascanfiles;     
end

%% specify which files to load for classification
if S.xval
    S.filenames = S.filenames_train;
else
%     S.filenames_h{1} = S.filenames_train;
%     S.filenames_h{2} = S.filenames_test;
    S.filenames = [S.filenames_train; S.filenames_test];
end
    S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);

%% Runs Parameters

%S.runs_vector - number of volumes per each run
if S.xval
    S.runs_vector =  [par.(par.task).numvols(S.TrainRuns)];
else
    S.runs_vector =  [parTr.(parTr.task).numvols(S.TrainRuns) par.(par.task).numvols(S.TestRuns)];
    S.TrainTestOrder = [ones(size(S.TrainRuns)) 2*ones(size(S.TestRuns))];
end


S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = par.TR;

%% Volume Parameters
S.vol_info = spm_vol(fullfile(par.subdir, 'TestExp/TestExp01/', ['swa', par.str, '_test01_0004.nii'])); %get functional data resolution info for spm .img writing

% % S.roiWithNonTaksVoxels = fullfile(par.anatdir, 'tnativeOccTemp.nii');
% % S.roiWithNonTaksVoxelsName = 'tnativeOccTemp.nii';

%S.roi_name = ['parietalNoPostCentral.img'];
%S.roi_file = '/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/parietalNoPostCentral.img';
S.roi_file = [S.expt_dir '/Masks/inv_SEPT09_MVPA_MASK_resliced4mm.nii'];
S.roi_name = 'inv_SEPT09_MVPA_MASK_resliced4mm.nii';

%S.roi_name = 'rc1V001.nii';
%S.roi_file = fullfile(S.anat_dir, S.roi_name);

%S.roi_name = 'tnativeOccTemp.nii';
%S.roi_file = fullfile(par.anatdir, S.roi_name);

% % S.noiseVoxels_file = fullfile(par.anatdir, 'rnativec2V001.nii');
% % S.noiseVoxels_name = 'rnativec2V001.nii';
% % 
% % S.sigVoxels_file = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
% % S.sigVoxels_name = 'tnativeOccTempGrey.nii';

S.secondaryMask = []; % secondary mask (the specific classification mask)
%masks the primary data loaded in the workspace. [] = no secondary mask.
%S.secondaryMask = fullfile(par.anatdir, 'tnativeOccTempGrey.nii');
%S.secondaryMask = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/OTnoHipp.img'];
%S.secondaryMask = ['/Users/gordonam/Studies/AG1/mvpa/OT_2012/rhipp.img'];
%S.secondaryMask = ['/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/masks/rparietal.img'];


%% Workspace Parameters
S.use_premade_workspace = 1;
%S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);



%% Pattern names
S.patternType = 'raw'; %'raw' or 'betas'
% S.preprocPatName = 'spiral_dn';
% S.preprocPatName = 'patsAllVox_z_dn';
S.preprocPatName = 'epi_d_hp_z';
S.preprocPatCondensedName = [S.preprocPatName '_condensed'];

if isempty(S.secondaryMask)
    S.preprocPatNameFinalMask = S.preprocPatName;
else
    S.preprocPatNameFinalMask = [S.preprocPatName '_masked'];
end

%% Artifacts
% % S.artFile = (fullfile(par.artrepdir, ['art_global_modified_' par.substr])); %directory where artifact information is contained
S.inactivateArtifacts = 0; %remove artifact trials?

%% Iteration Parameters
S.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data

%% Balancing Parameters
S.equate_number_of_trials_in_groups = 1; % equate number of trials in conditions 1 and 2
S.numBalancedIts = 1; % number of iterations to run, with different randomization for the balancing

%% Z-Scoring and outlier detection
S.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification 
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% Importance Maps
S.generate_importance_maps = 0; %visualize classifier weights
S.generateBetaMaps = 1; %use betas, instead of importance values
S.impType = {'pos' 'neg' 'both' 'raw'}; %importance map types
S.regNames = S.condsTrain;

%% Special types of analysis
S.searchlightAnalysis = 0; % run a searchlight analysis
S.linReg = 0; % run an analysis with a continuous outcome variable

%% Subsample
S.subsampleToMatch = 0; %subsample trials to match quantities across them.
S.balanceHiAndLowConf = 0;% match N of hi and low confidence trials?

%% voxel interactions
S.includeVoxelInteractions = 0; %include interactions among voxels?   
S.interactionType = 2;
S.intEst = 1;
S.intConcat = 0;
S.intPThresh = .001;
S.intReportIncrement = 100;
S.intFilePrefix = 'intVox';
S.intMaskName = 'interactionMask';
S.intPatName = 'interactions';
S.intGroupName = 'interactionsGroup';
S.intUseIntsWithUniqueInfo = 1;

%% Signal intensity analysis
S.thisSigIntenseSelector = 'randomNFold_xval'; %which selector to use for signal intensity analysis
S.zscoreIntensityVals = 1; % zscore the intensity values?

%% Denoising
S.denoise = 0; %undergo denoising?
S.denoiseOpt.denoisespec = '10001'; %which parts of the glm output do we want to save?

%% Mean Signal Extraction Params
% parameters for selecting the mean signal from a class-specific ROI for each pattern.

S.extractMeanSignal = 0; %1 - do signal extraction. 0 - don't do this. 
S.defineROIsFromANOVAFS = 0; % define ROIs using ANOVA-based feature selection, instead of pre-defining them. 
S.logreg_2Features = 0; %perform a logistic regression, using the two extracted intensity vectors

% % S.ROI1PatName = [S.preprocPatCondensedName '_ROI1'];
% % S.ROI1_name = ['occipitoTemporal_faceVsScene_500vox.img'];
% % S.ROI1_file  = [par.subdir '/analysis_loc_mnem/' S.ROI1_name];
% % 
% % S.ROI2PatName = [S.preprocPatCondensedName '_ROI2'];
% % S.ROI2_name  = ['occipitoTemporal_sceneVsFace_500vox.img'];
% % S.ROI2_file   = [par.subdir '/analysis_loc_mnem/' S.ROI2_name];

%% TR Weighting
%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier?  S.TR_weights_set{1} -
%training weights, S.TR_weights_set{2} = testing weights
S.TR_weights_set = {[0 0 0 .5 .5] [0 0 0 .5 .5]}; 

%% classifier parameters
S.class_args.train_funct_name = 'train_pLR'; %training function
S.class_args.test_funct_name = 'test_pLR'; %testing function
S.class_args.classType = 'pLR';
S.perfmet_functs = 'perfmet_maxclass'; % performance metric
S.statmap_funct = 'AG_statmap_anova'; 
S.nPlsCompsSet = 0; % number of pls components to include. 0 = do not use pls components.
S.nFolds = 10; % number of cross validation iterations

S.class_args.nVox = 0; % number of voxels to select with feature selection e.g. [1000 5000 10000]
S.class_args.libLin = '-q -s 0 -B 1'; %arguments for liblinear
S.class_args.libsvm = '-q -s 0 -t 2 -d 3'; % arguments for libsvm
S.class_args.constant = true; % include a constant term?
S.class_args.prefitWeights = true; 
S.class_args.chooseOptimalPenalty = 1; % cycle through cost parameters in the training set, and chose the optimal one?
S.class_args.penaltyRange = [.001 .005 .01 .05 .1 .5 1 5 10 50 100 500 1000 50000]; % cost parameters to cycle through
S.class_args.radialBasisSelection = [];%[.00001 .0001 .001 .01 .1 1 10];
S.class_args.nFoldsPenaltySelection = 10; % number of cross validation folds for penalty parameter selection. 

S.class_args.penalty = 1;

end
