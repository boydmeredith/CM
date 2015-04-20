function meanSmooth = cm_smoothness_for_searchlight(smoothest_imgtype, run_3dclustsim, tr_to_test)
if nargin < 1
    smoothest_imgtype = 'null_srch';
end
if nargin < 2
    run_3dclustsim = 1;
end

groupMask = '/Volumes/awagner/wagner/jtboyd/CMret/Data/Functional/Masks/SEPT09_MVPA_MASK_resliced4mm.nii';
mvpa_root_dir = '/Volumes/awagner/wagner/jtboyd/CMret/Data/Functional/sl_mvpa';
glm_root_dir = '/Volumes/awagner/wagner/thackery/CM/Data/' ;
s_array = [1 3:10 12:26];
smoothVec = nan(length(s_array), 3);
subj_eg_file = '/Volumes/awagner/wagner/jtboyd/CMret/Data/Functional/CM001/mvpa/CM001_SEPT09_MVPA_MASK_resliced4mm.nii_s8mm_wa.mat';
subj_eg = load(subj_eg_file);
subj_eg.header.id='group';
idx=1;

%iterate through subjects getting FWHM from the residual time series
for s = s_array
    id = sprintf('CM%03d',s);
    switch smoothest_imgtype
        case 'res_ms'
            test_dir = fullfile(glm_root_dir, id, 'Results_model2', 'Test');
            %get the square root of the mean squared residual
            system(['3dcalc -prefix ' fullfile(test_dir, 'sqrt_ResMS') ...
                ' -a ' fullfile(test_dir, 'ResMS.hdr') ' -expr ''sqrt(a)''']); %Convert ResMS.hdr file to an image more similar to the residuals images
            img_for_smoothest = fullfile(test_dir, 'sqrt_ResMS+tlrc');
            mask_for_smoothest = fullfile(test_dir, 'mask.hdr');
        case 'null_srch'
%             %if the user supplies multiple trs, recursively run the
%             %function, but skip the clustsim step at the end
%             if length(tr_to_test) > 1
%                 for trnum = 1:length(tr_to_test)
%                     meanSmooth(trnum,:) = cm_smoothness_for_searchlight(smoothest_imgtype, 0, tr_to_test(trnum));
%                 end
%                 return
%             end
            
            test_dir = fullfile(mvpa_root_dir, sprintf('trByTr_exToCm_hitVCrs_wholebrain_crossTr_tr%iToTr%i',tr_to_test, tr_to_test));
             %load this subject's searchlight results file
            sr = load(fullfile(test_dir, [id '.mat']));
            %write the subjects scrambled regs to an spm image
            subj = subj_eg.subj;
            subj = initset_object(subj, 'pattern', ...
                ['null_srch_' id], ...
                sr.res.auc_scram_regs,'masked_by','SEPT09_MVPA_MASK_resliced4mm.nii')
            subj.patterns{end}.header.vol = subj.patterns{end-1}.header.vol{1}
            subj.patterns{end}.header.vol.fname = fullfile(test_dir,['null_srch_' id '.img']);
            write_to_spm(subj,'pattern',['null_srch_' id ]);
            img_for_smoothest = fullfile(test_dir, ['null_srch_' id '.hdr']); %use the image of aucs under the null hypothesis
            mask_for_smoothest = groupMask; %use the mask applied to the searchlight image
        otherwise
            error('Error: can''t specified type of image for smoothness calculation\n');
    end
    smooth_fname = fullfile(test_dir, [id '_smoothness.txt']);

    %get the FWHM of residual and stick it in <smoothness.txt>
    system(['3dFWHMx -arith -mask ' mask_for_smoothest ...
        ' -dset ' img_for_smoothest ...
        ' -out ' smooth_fname]);
    
    %this is what it would like if we did this in spm, but it doesn't
    %like the resulting file formats so doing it afni is fine
    %smoothVec(idx, :) = spm_est_smoothness(fullfile(test_dir, 'sqrt_ResMS+tlrc.BRIK'), fullfile(test_dir, 'mask.hdr'));
    
    %read the values in from the output text file
    fid = fopen(smooth_fname);
    smoothVec(idx,:) = cell2mat(textscan(fid, '%f %f %f'));
    
    %increment
    idx = idx+1;
end
%get the mean of the fwhm across subjects
meanSmooth = mean(smoothVec);

if run_3dclustsim
    system(['3dClustSim -pthr .05 .02 .01 .005 .001  -mask ' groupMask ' -fwhmxyz ' num2str(meanSmooth(1)) ' ' num2str(meanSmooth(2)) ' ' num2str(meanSmooth(3))]);
end