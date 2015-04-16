root_dir = '/Volumes/awagner/wagner/jtboyd/CMret/Data/Functional/sl_mvpa';
cd(root_dir);
tr_sl_dirs = dir('trByTr*');
tr_sl_dir_names = {tr_sl_dirs.name};
file_to_test = 'abs_mean_auc_diff_chance_p0.01_minuspt5.nii'
extents    = [53 63 86]; % minimal cluster size.. 
%iterate through result images
for d = 1:length(tr_sl_dir_names)
    %iterate through cluster thresholds
    for k = extents
        cd(fullfile(root_dir, tr_sl_dir_names{d}));
        fprintf('Examining results from %s\nApplying extent threshold of %d\n', tr_sl_dir_names{d}, k);
        ROI = fullfile(root_dir, tr_sl_dir_names{d}, file_to_test); % input image (binary, ie a mask)
        ROIf = fullfile(root_dir, tr_sl_dir_names{d}, sprintf('%s_k%04d.nii',file_to_test,k)); % output image (filtered on cluster size)
        %-Connected Component labelling
        V = spm_vol(ROI);
        dat = spm_read_vols(V);
        [l2, num] = spm_bwlabel(double(dat>0),26);
        if ~num, warning('No clusters found.'); end
        %-Extent threshold, and sort clusters according to their extent
        [n, ni] = sort(histc(l2(:),0:num), 1, 'descend');
        l  = zeros(size(l2));
        n  = n(2:end); ni = ni(2:end)-1;
        ni = ni(n>=k); n  = n(n>=k);
        for i=1:length(n), l(l2==ni(i)) = i; end
        clear l2 ni
        fprintf('Selected %d clusters (out of %d) in image.\n',length(n),num);

        %-Write new image
        V.fname = ROIf;
        spm_write_vol(V,l~=0); % save as binary image. Remove '~=0' so as to
                           % have cluster labels as their size. 
                           % or use (l~=0).*dat if input image was not binary
    end
end


 