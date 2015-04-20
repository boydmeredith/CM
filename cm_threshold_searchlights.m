root_dir = '/Volumes/awagner/wagner/jtboyd/CMret/Data/Functional/sl_mvpa';
cd(root_dir);
tr_sl_dirs = dir('trByTr*');
tr_sl_dir_names = {tr_sl_dirs.name};
file_to_test = 'mean_auc_diff_chance_p0.01_minuspt5.nii'
extents    = [22 29]; % minimal cluster size..
%iterate through result images
for d = 1:length(tr_sl_dir_names)
    %iterate through cluster thresholds
    for k = extents
        cd(fullfile(root_dir, tr_sl_dir_names{d}));
        fprintf('Examining results from %s\nApplying extent threshold of %d\n', tr_sl_dir_names{d}, k);
        ROI = fullfile(root_dir, tr_sl_dir_names{d}, file_to_test); % input image (binary, ie a mask)
        ROIf = fullfile(root_dir, tr_sl_dir_names{d}, sprintf('both_%s_k%04d.nii',file_to_test,k)); % output image (filtered on cluster size)
%         ROIf_pos = fullfile(root_dir, tr_sl_dir_names{d}, sprintf('pos_%s_k%04d.nii',file_to_test,k)); % output image (filtered on cluster size)
%         ROIf_neg = fullfile(root_dir, tr_sl_dir_names{d}, sprintf('neg_%s_k%04d.nii',file_to_test,k)); % output image (filtered on cluster size)
%         %-Connected Component labelling
%         V = spm_vol(ROI);
%         dat = spm_read_vols(V);
%         [l2_pos, num_pos] = spm_bwlabel(double(dat > 0),26);
%         
%         %% positive values:
%         if ~(num_pos), warning('No positive clusters found.'); end
%         %-Extent threshold, and sort clusters according to their extent for
%         %positive values
%         [n_pos, ni_pos] = sort(histc(l2_pos(:),0:num_pos), 1, 'descend');
%         l_pos  = zeros(size(l2_pos));
%         n_pos  = n_pos(2:end); ni_pos = ni_pos(2:end)-1;
%         ni_pos = ni_pos(n_pos>=k); n_pos  = n_pos(n_pos>=k);
%         for i=1:length(n_pos), l_pos(l2_pos==ni_pos(i)) = i; end
%         clear l2_pos ni
%         fprintf('Selected %d positive clusters (out of %d) in image.\n',length(n_pos),num_pos);
%         %-Write new image
%         V.fname = ROIf_pos;
%         spm_write_vol(V,l_pos~=0); % save as binary image. Remove '~=0' so as to
%         % have cluster labels as their size.
%         % or use (l~=0).*dat if input image was not binary
%         
%         %% negative values:
%         [l2_neg, num_neg] = spm_bwlabel(double(dat < 0),26);
%         if ~(num_neg), warning('No negative clusters found.'); end
%         %-Extent threshold, and sort clusters according to their extent for
%         %negative values
%         [n_neg, ni_neg] = sort(histc(l2_neg(:),0:num_neg), 1, 'descend');
%         l_neg  = zeros(size(l2_neg));
%         n_neg  = n_neg(2:end); ni_neg = ni_neg(2:end)-1;
%         ni_neg = ni_neg(n_neg>=k); n_neg  = n_neg(n_neg>=k);
%         for i=1:length(n_neg), l_neg(l2_neg==ni_neg(i)) = i; end
%         clear l2_neg ni
%         fprintf('Selected %d negative clusters (out of %d) in image.\n',length(n_neg),num_neg);
%         %-Write new image
%         V.fname = ROIf_neg;
%         spm_write_vol(V,l_neg~=0); % save as binary image. Remove '~=0' so as to
%         % have cluster labels as their size.
%         % or use (l~=0).*dat if input image was not binary
        
        
        
        %% both positive and negative values:
        [l2, num] = spm_bwlabel(double(dat < 0 | dat > 0),26);
        if ~(num), warning('No negative or positive clusters found.'); end
        %-Extent threshold, and sort clusters according to their extent for
        %negative values
        [n, ni] = sort(histc(l2(:),0:num), 1, 'descend');
        l  = zeros(size(l2));
        n  = n(2:end); ni = ni(2:end)-1;
        ni = ni(n>=k); n  = n(n>=k);
        for i=1:length(n), l(l2==ni(i)) = i; end
        clear l2 ni
        fprintf('Selected %d negative and positive clusters (out of %d) in image.\n',length(n),num);
        %-Write new image
        V.fname = ROIf;
        spm_write_vol(V,(l~=0).*dat); % save as binary image. Remove '~=0' so as to
        % have cluster labels as their size.
        % or use (l~=0).*dat if input image was not binary
    end
end


