function subjects_tc = CM_batch_FIR_marsbar(contrast, saveit)


%calculates FIR for specified ROI for all subjects
%MUST HAVE DEFINED ROI IN MARSBAR PRIOR TO RUNNING THIS SCRIPT
%(open marsbar, 'ROI definition', call 2nd level SPM.mat of interest,'Write ROIs')
%THEN MODIFY 'ROI_FILE' VARIABLE IN THIS SCRIPT AS NEEDED FOR EACH ROI
%written by Melina Uncapher 5/23/08, adapted from Matthew Brett's scripts
if nargin < 1
    contrast = 'Task Effects (Explicit Memory Runs_CM Runs) x (Hits_CRs)';
end

if nargin < 2
    saveit = 1;
end

FIR_LENGTH = 16;

FMRI_DIR = '/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional';
ROI_DIR = fullfile(FMRI_DIR, 'Masks');
GROUP_DIR = '/Users/Jesse/fMRI/COUNTERMEASURES/Group_analyses_univariate';
MODEL_DIR = fullfile(GROUP_DIR, 'Model2/Test');
SAVE_NAME = 'subjects_tc_rois.mat';
CONDS = {'AA_Explicit_Hit' 'EA_Explicit_Hit' 'AA_CM_Hit' 'EA_CM_Hit'...
    'AA_Explicit_CR' 'EA_Explicit_CR' 'AA_CM_CR' 'EA_CM_CR'};
SUBJ_ARR = [1,3:10,12:26];
subjno_to_id = @(x) sprintf('CM%03d',x);
subjno_to_analysisdir = @(x) fullfile(FMRI_DIR, subjno_to_id(x), 'Results_model2/Test');

analysis_dir = fullfile(MODEL_DIR, contrast);

SUBJ_ARR = 1


roi_names = {'rleftMtlMask_tbm112514.nii'}%, 'rleftSplMask_tbm112514.nii', 'rleftAngMask_tbm112514.nii'};

for isubj=1:length(SUBJ_ARR)
    subjno = SUBJ_ARR(isubj);
    
    sub_dir = subjno_to_analysisdir(subjno);
    
    spm_name = fullfile(sub_dir, 'SPM.mat');
    
    D  = mardo(spm_name);
    
    D = autocorr(D,'fmristat',2);
    
    for jclus = 1:length(roi_names)
        
        % Make marsbar design object
        roi_file = fullfile(ROI_DIR, roi_names{jclus});
        
        % Make marsbar ROI object
        R  = maroi_image(roi_file);
        
        % Fetch data into marsbar data object
        Y  = get_marsy(R, D, 'mean');
        
        % Get contrasts from original design
        xCon = get_contrasts(D);
        
        % Estimate design on ROI data
        E = estimate(D, Y);
        
        %---------------------
        %NOW FIR...
        % Get definitions of all events in model
        [e_specs, e_names] = event_specs(E);
        n_events = size(e_specs, 2);
        
        % Bin size in seconds for FIR
        bin_sz = tr(E);
        
        % Length of FIR in seconds
        FIR_LENGTH = 16;
        
        % Number of FIR time bins to cover length of FIR
        nbins = FIR_LENGTH / bin_sz;
        
        % Options - here 'single' FIR model, return estimated % signal change
        opts = struct('single', 1, 'percent', 1);
        
        fir_tc = [];
        % Return time courses for all events in fir_tc matrix
        for e_s = 1:n_events
            fir_tc(:, e_s) = event_fitted_fir(E, e_specs(:,e_s), bin_sz, ...
                nbins, opts);
        end
        
        for kcond = 1:length(CONDS)
            e_s = find(strcmp(CONDS{kcond}, e_names));
            
            if isempty(e_s)
                fir_tc(:, kcond) = nan(fir_length/bin_sz,1);
                subjects_tc(isubj).rois(curr_clust).maxWithNeighbors(kcond) = NaN;
            else
                
                fir_tc(:, kcond) = event_fitted_fir(E, e_specs(:,e_s), bin_sz, nbins, opts);
                
                [maxC, maxI] = max(fir_tc(2:((FIR_LENGTH/bin_sz)-1),kcond));
                
                
                %subjects_tc(isubj).rois(jclus).maxval(c) = maxC;
                
                subjects_tc(isubj).rois(jclus).maxval(kcond) = maxC;
                neighbors = intersect(1:size(fir_tc,2), [maxI-1: maxI+1]);
                subjects_tc(isubj).rois(jclus).maxWithNeighbors(kcond) = mean(fir_tc(neighbors,kcond));
                
            end
            subjects_tc(isubj).rois(jclus).conds = CONDS;
            
            subjects_tc(isubj).rois(jclus).roiName = roi_names{jclus};
            subjects_tc(isubj).rois(jclus).raw = fir_tc;
        end
        
        save(fullfile(analysis_dir, SAVE_NAME), 'subjects_tc');
        
        
    end
    subjects_tc(isubj).conds = CONDS;
    subjects_tc(isubj).e_specs = e_specs;
    subjects_tc(isubj).e_names = e_names;
    subjects_tc(isubj).n_events = n_events;
    
end






