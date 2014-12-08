function betas_out = CM_extract_betas(clus_coords, csvName)
% in the intended use case <clus_coords> is created using spm's eigenvariate
% function, which outputs a variable <xY.XYZmm>

USING_TD_DERIVS=1;
% experiment specific parameters
expt_name = 'CM';
fmri_dir = '/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional';
analysis_name = 'Results_model2';
saveDir = '/Users/Jesse/fMRI/COUNTERMEASURES/Group_analyses_univariate/Model2/Test/Task Effects (Explicit Memory Runs_CM Runs) x (Hits_CRs)';
subj_array = [1,3:10,12:26];

display('Warning: this function will not correct for instances where subjects do not share conditions');

% get count of subjects
nSubj = length(subj_array); 

% ensure that clus_coords is properly oriented
if size(clus_coords,2) ~= 3
    display('reshaping input voxel coordinates');
    clus_coords = clus_coords';
end
[nRows_vox, nCols_vox] = size(clus_coords);
assert(nCols_vox == 3);


for vox = 1:nRows_vox
    
    fprintf('\n Pulling betas from voxel %d of %d ',...
               vox, nRows_vox);
    
    this_vox = clus_coords(vox,:);
    
    for subj = 1:nSubj
        
        this_subjname = sprintf('%s%03d', expt_name, subj_array(subj));
               
        this_subj_analysis_dir = fullfile(fmri_dir, this_subjname, analysis_name, 'Test');
        
        this_subj_conds = load(fullfile(this_subj_analysis_dir,  sprintf('CM%03d_Test_Onsets.mat',subj_array(subj))), 'names');
        
        nConds = length(this_subj_conds.names);
        
        for cond = 1:nConds
            
            this_condname = this_subj_conds.names{cond};
            
            %choose every third beta, assuming that we're using time and
            %dispersion derivatives
            if USING_TD_DERIVS
                this_beta_fname = fullfile(this_subj_analysis_dir, sprintf('beta_%04d.img',1+(cond-1)*3));
            else
                this_beta_fname = fullfile(this_subj_analysis_dir, sprintf('beta_%04d.img',cond));
            end
            this_beta = get_vox(this_beta_fname, this_vox);
            
            betas_out.(this_condname)(subj,vox) = this_beta;
            
            if isnan(this_beta)
                error('beta was returned as nan');
            end
            

        end
        
    end
    
end


    fn = fieldnames(betas_out);
    for f = 1:length(fn)
        out.(fn{f}) = vertcat(mean(betas_out.(fn{f}),2));
        toPrint(:,f) = [fn{f}; num2cell(out.(fn{f}))];
    end
    
    if nargin>1
        cell2csv(fullfile(saveDir, [csvName '.csv']), toPrint, ',', 2000);
    end
end
            

