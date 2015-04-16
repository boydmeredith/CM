function [auc_80_bins, fas_80, hits_80] = getAuc(out,ix,plotit)
% sort by signed classifier "confidence" (for ROI curves)
actsVec = out.actsVec(ix);
correctsVec = out.correctsVec(ix);
desiredsVec = out.desiredsVec(ix);
[sorted_diffs ind] = sort(actsVec,'descend');
correct_sorted = correctsVec(ind);
desireds_sorted = desiredsVec(ind);

for i = 1:length(sorted_diffs);
    hit_rate(i) = length(intersect(find(desireds_sorted == 1),[1:i])) / length(find(desireds_sorted == 1));
    fa_rate(i) = length(intersect(find(desireds_sorted == 2),[1:i])) / length(find(desireds_sorted == 2));
end
auc = auroc(hit_rate',fa_rate');

% create ROC function with 80 discrete bins (this allows the ROC curves to be more easily averaged across individuals)
roc_bin_intervals = .975:-.025:-1;
for bin_num = 1:80
    hits_80(bin_num)=length(intersect(find(desireds_sorted == 1),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 1));
    fas_80(bin_num)=length(intersect(find(desireds_sorted == 2),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 2));
end
auc_80_bins = auroc(hits_80',fas_80');

if plotit    
    plot(fa_rate,hit_rate,'.-');
    hold on;
    plot([0 1],[0 1],'r');
    set(gcf,'Color','w');
    %     xlabel('P(Old|New)')
    %     ylabel('P(Old|Old)')
end
end
