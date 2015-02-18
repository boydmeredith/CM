function [expt, classArgs, par]= CM_mvpa_params(subj_id, task, proclus)

% establish parameters for mvpa analysis
% <subj_id> - identifier for the given subject. Can be numerical or a
% string
% <task> 'perc' or 'mnem'

%% EXPT SPECIFIC INFO
expt.name = 'CM_ret';
if proclus
	expt.dir = '/hsgs/projects/awagner/jtboyd/CM_ret/Data/Functional/';
else
    expt.dir = '/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional/'; % top-level expt. directory where individual subject folders live
end
expt.mvpaDirStr = 'mvpa';
%expt.dataImgsToUse = 'raw_filenames.mat'; % .mat file containing names of all functional images (must exist for each subject; can be created by running cellstr(SPM.xY.P) on subject's SPM.mat file)
expt.numTpPerRun = 256; % number of TRs per scanning run (coded for fixed value; adjust code structure if variable TR counts across runs)
%expt.roiName = 'rLR_PrePostCentralandSMA';  % name of mask to be used for voxel selection (can be a small ROI, a whole-brain mask, or anywhere in between)
expt.xval_iters_for_imp_map = [];
%flags determine how classifier wrapper behaves
expt.num_full_iter = 1; % number of times to run the entire classification process, including feature selection
expt.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
expt.num_results_iter = 1; % number of times to run the post-feature selection classification process (select subset of the data and train/test classifier)
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
expt.roiName = 'SEPT09_MVPA_MASK_resliced4mm.nii';
%expt.roiName = 'rleftSplMask_tbm112514.nii';
expt.trnSubCondsToBal = {};
expt.tstSubCondsToBal = {};

expt.trsToAverageOver = [1 2 3 4 5 6];
expt.trWeights_train = [0 0 .33 .34 .33 0];
expt.trWeights_test = [0 0 .33 .34 .33 0];
expt.nVox_thresh = 0;
expt.subjIdFormat = 'CM%03d';
expt.onsetsFname = 'mvpa_ons.mat';
expt.trSecs = 2; %used to convert from seconds to TRs
expt.perfmetFuncts = 'perfmet_maxclass';
expt.group_mvpa_dir = [expt.dir '/mvpa_results'];
expt.which_traintest = [1 2 3 4 0 0 0 0];
expt.condNames = {'hit', 'cr'};
expt.subjFname = fullfile(expt.dir, subj_id, expt.mvpaDirStr,[subj_id '_' expt.roiName '.mat']); %having this previously saved avoids time-consuming data extraction and preproc

expt.srch_radius = sqrt(2);


%% establish general parameters

par = CM_Params(subj_id, task, 1, proclus);

end
