function auc  = calculate_auc(acts_vector, desireds_vector, plotit)
% auc  = calculate_auc(acts_vector, desireds_vector)

if nargin < 3
plotit =0
end


acts_class1 = acts_vector(1,:);
[sorted_diffs ind] = sort(acts_class1,2,'descend')
desireds_sorted = desireds_vector(ind);
for i = 1:length(sorted_diffs)	
hit_rate(i) = length(intersect(find(desireds_sorted==1),[1:i])) / length(find(desireds_sorted == 1));
fa_rate(i) = length(intersect(find(desireds_sorted==2),[1:i])) / length(find(desireds_sorted == 2));
end

auc = auroc(hit_rate', fa_rate');

if plotit
plot(fa_rate, hit_rate)
end
