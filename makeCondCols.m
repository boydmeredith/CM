function condCols = makeCondCols(expt, names)
%% sets the indices of the conditions that need to be examined
condCols = zeros(1,length(expt.condNames));
for i = 1:length(condCols)
    thisName = expt.condNames{i};
    condCols(i) = find(ismember(names,thisName));
end
end
