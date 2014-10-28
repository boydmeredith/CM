function res = CM_run_mvpa_srch(subjArray, varargin)


% this selector will be set up by cm_create_sl_selectors and its
% name will be passed to feature_select
SL_SELECTOR_TO_USE = 'trEx_teCm_sl';

for iSub = subjArray
	[Expt, classArgs, ~, par] = CM_mvpa_params(iSub, 'ret');
	Expt.sl_selector_to_use = SL_SELECTOR_TO_USE;
	subjId = sprintf(Expt.subjIdFormat,iSub);
	Expt.saveName = cm_make_sl_savename(Expt);
	Expt.impMapDirStr = fullfile(Expt.dir, '/mvpa_results/importance_maps/', Expt.saveName);

	Expt.roiFname = fullfile(Expt.dir,'Masks',Expt.roiName);
	thisMvpaDir = fullfile(Expt.dir, subjId, Expt.mvpaDirStr);
	Expt.subjFname = fullfile(thisMvpaDir, [subjId '_' Expt.roiName '_s8mm_wa.mat']);
	if ~exist(Expt.subjFname,'file')
		Subj = CM_mvpa_load_and_preprocess_raw_data(subjId, Expt, nRuns, 1);
	else 
		Subj = load(Expt.subjFname,'subj');
		Subj=Subj.subj;
	end

	onsetsFile = fullfile(thisMvpaDir, Expt.onsetsFname);
	Expt.ons = load(onsetsFile);
	Expt.scanfiles = vertcat(par.swascanfiles.(par.task));

	fprintf('\n\nPrepping onsets and patterns for CM%03d...',iSub);
	Subj = cm_condense_onsets_and_patterns(Subj, Expt);

	fprintf('\n\nPrepping selectors for CM%03d...',iSub);
	Subj = cm_create_sl_selectors(Subj);
	for iResIteration = 1:Expt.num_results_iter
		new_active_trials = [];
		if Expt.equate_number_of_trials_in_cond_1_and_2 
			Subj = create_balanced_xvalid_selectors(Subj, 'conds',SL_SELECTOR_TO_USE);
		end
		nSelectors = length(find_group(Subj,'selector',[SL_SELECTOR_TO_USE '_bal']));
		for iSelector = 1:nSelectors
			new_active_trials = horzcat(new_active_trials, find(Subj.selectors{end-nSelectors+iSelector}.mat==2));
		end
		new_actives_selector(new_active_trials)=1;
		Subj= initset_object(Subj,'selector','conditions_of_interest_bal_within_runs', new_actives_selector);
		active_trials = new_active_trials;

		if Expt.perform_second_round_of_zscoring
			display('Performing second round of z-scoring');
			Subj.patterns{5}.mat(:,active_trials) = zscore(Subj.patterns{5}.mat(:,active_trials)')';
		end

		srch_radius = Expt.srch_radius;
		Subj.adj_sphere = create_adj_list(Subj,Expt.roiName,'radius',srch_radius);

		class_args.train_funct_name = 'train_gnb';
		class_args.test_funct_name = 'test_gnb';

		scratch.class_args = class_args;
		scratch.perfmet_funct = 'perfmet_maxclass';
		scratch.perfmet_args = struct([]);

		statmap_srch_arg.adj_list = Subj.adj_sphere;
		statmap_srch_arg.obj_funct = 'statmap_classify';
		statmap_srch_arg.scratch = scratch;

		Subj = feature_select(Subj, ...
			'epi_d_hp_z_condensed', ... %data
			'conds', ... % binary regs (for gnb)
			SL_SELECTOR_TO_USE, ... % selector
			'statmap_funct', 'statmap_searchlight', ... %function
			'statmap_arg',statmap_srch_arg, ...
			'new_map_patname', 'epi_d_hp_z_condensed_srch', ...
			'thresh', []);

		
        
        
        Subj=create_sorted_mask( ...
            Subj,'epi_d_hp_z_condensed_srch', ...
            'epi_d_hp_z_condensed_srch_top200', ...
            200,'descending',true)
        
        Subj=create_sorted_mask( ...
            Subj,'epi_d_hp_z_condensed_srch', ...
            'epi_d_hp_z_condensed_srch_bottom200', ...
            200,'descending',false)
        
        class_args.train_funct_name = 'train_plr';
        class_args.test_funct_name = 'test_plr';
        class_args.penalty = 10;
        
        [Subj res.top200SlMask{iSub}] = cross_validation( ...
            Subj, ...
            'epi_d_hp_z_condensed', ... % data
            'conds', ... % regressors
            SL_SELECTOR_TO_USE, ... % xval selectors
            'epi_d_hp_z_condensed_srch_top200', ... % mask group
            class_args,'perfmet_functs', expt.perfmetFuncts);
        
        [Subj res.bottom200SlMask{iSub}] = cross_validation( ...
            Subj, ...
            'epi_d_hp_z_condensed', ... % data
            'conds', ... % regressors
            SL_SELECTOR_TO_USE, ... % xval selectors
            'epi_d_hp_z_condensed_srch_bottom200', ... % mask group
            class_args,'perfmet_functs', expt.perfmetFuncts);
        

	end

end

if ~exist(Expt.groupMvpaDir)
	mkdir(Expt.groupMvpaDir);
end
save(fullfile(Expt.group_mvpa_dir, Expt.saveName),'res');

end

function name = cm_make_sl_savename(Expt)
condNamesStr = strrep(strjoin(Expt.condNames),',','V');
trTeStr = strrep(num2str(Expt.which_traintest),'  ','_');
trWStr = strrep(num2str(Expt.trWeights),' ','_');
roiStr = regexprep(Expt.roiName,'(\w*).*','$1');
name = sprintf('srchlight_onds_%s_trTe%s_trW%s_roi%s',condNamesStr,trTeStr,trWStr,roiStr);

end

function Subj = cm_condense_onsets_and_patterns(Subj, Expt);   
condensed_runs = [];
nExamples = [];

names = Expt.ons.names;
onsets = Expt.ons.onsets;

nScanFiles = length(Expt.scanfiles);
nRuns = nScanFiles/Expt.numTpPerRun;
nPatts = size(Subj.patterns{end}.mat,2)
nCondsTotal = size(onsets,2);
assert(nScanFiles == nPatts);

Expt.condCols = makeCondCols(Expt, names);
nCondsToClassify = length(Expt.condNames);
condOnsets = zeros(nCondsToClassify, nPatts);
i = 1;
for iCond = Expt.condCols;
	nExamples(i) = length(onsets{iCond});
	time_idx = onsets{iCond}/Expt.trSecs + 1;
	condOnsets(i, time_idx) = 1;
	i = i + 1;	
end

restTp = (sum(condOnsets,1) == 0);
condensedCondRegs = condOnsets(:,~restTp);

all_trials = sum(condOnsets,1);
for trNum = Expt.trsToAverageOver
	data_by_tr(trNum,:,:) = Expt.trWeights(trNum)*Subj.patterns{end}.mat(:,find(all_trials)+trNum-1);
end
temporally_condensed_data = squeeze(sum(data_by_tr(Expt.trsToAverageOver,:,:),1));
clear data_by_tr;

if Expt.remove_outlier_trials
	mean_across_voxels = mean(temporally_condensed_data,1);
	z_mean_across_voxels = zscore(mean_across_voxels);
	upper_outliers = find(z_mean_across_voxels> Expt.remove_outlier_trials);
	lower_outliers = find(z_mean_across_voxels< -1 * Expt.remove_outlier_trials);
	all_outliers = union(upper_outliers,lower_outliers)
	condensedCondRegs(:,all_outliers) = 0;
end

condensed_runs = Subj.selectors{1}.mat(~restTp);

Subj = initset_object(Subj, 'regressors', 'conds', condensedCondRegs, 'condnames', Expt.condNames);
Subj = duplicate_object(Subj, 'pattern', 'epi_d_hp_z','epi_d_hp_z_condensed');
Subj = set_mat(Subj, 'pattern', 'epi_d_hp_z_condensed', temporally_condensed_data,'ignore_diff_size',true);
Subj = add_history(Subj,'pattern','epi_d_hp_z_condensed','Pattern created created by JR custom code')
Subj = remove_mat(Subj,'pattern','epi_d_hp_z');

Subj.selectors{1}.mat = condensed_runs;
Subj.selectors{1}.matsize = size(condensed_runs);

active_trials = find(sum(condensedCondRegs));
actives_selector = zeros(1,size(condensedCondRegs,2));
actives_selector(active_trials) =1;
Subj = initset_object(Subj,'selector','conditions_of_interest',actives_selector);

end


function Subj = cm_create_sl_selectors(Subj)
FIRST_CM_RUN = 5;
CM_RUNS = 5:8;
runs = get_mat(Subj,'selector','runs');
nRuns = max(runs);
nTimepoints = length(runs);

runs_xval_sl = ones(nRuns, nTimepoints);
% set up train explicit test countermeasures selector by setting all explicit trials to be training examples (1s), 3 of 4 CM trials are set as searchlight training runs (3s), and the remaining CM trial is set as a final generalization test trial (2)
trEx_teCm_sl = ones(4, nTimepoints);
cm_run_ix = find(ismember(runs,CM_RUNS));
trEx_teCm_sl(:, cm_run_ix) = 3;
for r = FIRST_CM_RUN:nRuns
    r_ix = r-FIRST_CM_RUN+1;
    cur_final_test_run = find(runs == nRuns-r_ix+1);
	trEx_teCm_sl(r_ix, cur_final_test_run) = 2;
	cur_name = sprintf('trEx_teCm_sl_%i',r_ix);
	Subj = initset_object(Subj, 'selector', cur_name, ...
		trEx_teCm_sl(r_ix,:),'group_name','trEx_teCm_sl');
end

imagesc(trEx_teCm_sl);
set(gca, 'CLim', [1 3]);
colorbar;

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
	cur_name = sprintf('runs_xval_sl_%i',r);
	Subj = initset_object(Subj, 'selector', cur_name, ...
		runs_xval_sl(r,:), ...
		'group_name', 'runs_xval_sl' );
end

%imagesc(runs_xval_sl)
%colorbar



end


