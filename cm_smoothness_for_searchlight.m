run_3dclustsim = 0
groupMask = '/Volumes/awagner/wagner/jtboyd/CMret/Data/Functional/Masks/SEPT09_MVPA_MASK_resliced4mm.nii';
root_dir = '/Volumes/awagner/wagner/thackery/CM/Data/';
s_array = [1,3:10,12:26]
smoothVec = nan(length(s_array), 3);
idx=1;

%iterate through subjects getting FWHM from the residual time series
for s = s_array
    id = sprintf('CM%03d',s);
    test_dir = fullfile(root_dir, id, 'Results_model2', 'Test');
    cd(test_dir);
    
    %get the square root of the mean squared residual
    system(['3dcalc -prefix sqrt_ResMS -a ResMS.hdr -expr ''sqrt(a)''']); %Convert ResMS.hdr file to an image more similar to the residuals images
    
    %get the FWHM of residual and stick it in <smoothness.txt>
    system(['3dFWHMx -arith -mask ' fullfile(test_dir, 'mask.hdr') ...
        ' -dset ' fullfile(test_dir, 'sqrt_ResMS+tlrc') ...
        ' -out ' fullfile(test_dir, 'smoothness.txt')]);
    
    %this is what it would like if we did this in spm, but it doesn't
    %like the resulting file formats so doing it afni is fine 
    %smoothVec(idx, :) = spm_est_smoothness(fullfile(test_dir, 'sqrt_ResMS+tlrc.BRIK'), fullfile(test_dir, 'mask.hdr'));
    
    %read the values in from the output text file
    fid = fopen(fullfile(test_dir, 'smoothness.txt'));
    smoothVec(idx,:) = cell2mat(textscan(fid, '%f %f %f'));
    
    %increment 
    idx = idx+1;
end
%get the mean of the fwhm across subjects
meanSmooth = mean(smoothVec);

if run_3dclustsim
    system(['3dClustSim -pthr .05 .02 .01 .005 .001  -mask ' groupMask ' -fwhmxyz ' num2str(meanSmooth(1)) ' ' num2str(meanSmooth(2)) ' ' num2str(meanSmooth(3))]);
end