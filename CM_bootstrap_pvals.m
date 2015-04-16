function [aucdist pdist A B] = CM_bootstrap_pvals(nres_it, n_iter_samp, n_boot_samp, res, xval_it_toreport)


nsub = length(res.subjArray);
A = nan(nsub, n_boot_samp);


for isub = 1:nsub
    for jsamp = 1:n_boot_samp
        subNo = res.subjArray(isub);
        sub_res = res.subj{subNo};
        iter_to_sample = randi(nres_it,1,n_iter_samp);
        B(isub,jsamp,:) = iter_to_sample;
        for kiter = 1:length(iter_to_sample)
            sub_res_iter = sub_res.penalty.nVox.weights.iter{iter_to_sample(kiter)};
            A(isub,jsamp,kiter) = compute_auc(sub_res_iter,xval_it_toreport);
        
        end
    end
end

meanA = mean(A,3);
aucdist = mean(meanA);
[~, pdist] = ttest(meanA-.5);
hist(pdist)
hold on
line([.05/6 .05/6], [0 300])

    


function auc = compute_auc(sub_iter_res,xval_it_toreport) 
desireds_vector = []; prob_class_A_vector = []; 
for a = xval_it_toreport % concatenate the results across all cross-validation testing sets
    desireds_vector = horzcat(desireds_vector,sub_iter_res.iterations(a).perfmet.desireds); % create a vector of the true class of each trial (i.e., the desired outcome) (1=Class A; 2=Class B)
    prob_class_A_vector = horzcat(prob_class_A_vector, sub_iter_res.iterations(a).acts(1,:)); % create a vector of the classifier's scalar probability estimates that each trial is a trial of Class A
end

% sort by signed classifier "confidence" (for ROI curves)
[sorted_diffs ind] = sort(prob_class_A_vector,2,'descend');
desireds_sorted = desireds_vector(ind);

% create continuous ROC function
for i = 1:length(sorted_diffs);
    hit_rate(i) = length(intersect(find(desireds_sorted == 1),[1:i])) / length(find(desireds_sorted == 1));
    fa_rate(i) = length(intersect(find(desireds_sorted == 2),[1:i])) / length(find(desireds_sorted == 2));
end

auc = auroc(hit_rate',fa_rate')



end
end