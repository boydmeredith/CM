function [resB cvB out] = CM_batchMvpaPostProc(resDir, classStr, xvalIterToReport, task, subjArray, saveName) 
%function [resB cvB] = batchMvpaPostProc(dir, classStr, xvalIterToReport, task) 
%inputs:
%<resDir> the directory containing the results mat files
%<classSt>: default value is '*conds*". Filters mat files
%<xvalIterToReport>: allows specification of an array of xval iterations to report. This is specifically for cases when the xval iterations are very different, such as when one xval bin is all explicit trials and one xval bin is all countermeasures trials
%<task>: default value is 'ret'

toR = [];

if isempty(subjArray)
    subjArray = [1 3:10 12:26];
end

if isempty(classStr)
	classStr = '*conds*';
end
if isempty(task)
	task = 'ret';
end

dFN = dir(fullfile(resDir, [classStr '.mat']));

fileNames = {dFN.name};
dotFiles = cell2mat(cellfun(@(x) x(1)=='.', fileNames, 'UniformOutput', false));
fileNames(find(dotFiles)) = []; %remove hidden files that are prepended with dots.

for resf =1:length(fileNames)
    load(fileNames{resf});
    if length(res.subjArray) > length(res.subj)
        display(['skipping over ' res.subj{1}.penalty.nVox.weights.expt{1}.saveName ' b/c the results struct is shorter than the subject array...\n']);
        continue
    end
    [resB{resf} cvB{resf}] = CM_mvpaPostProc(res, xvalIterToReport, task, 1, fullfile(resDir,'figs'), subjArray);
    resB{resf}.name = fileNames{resf};
    cvB{resf}.name = fileNames{resf};
    
    toR.resname{resf} = resB{resf}.name;
    toR.mean_auc(resf,:) = nanmean(resB{resf}.auc);
    toR.sem_auc(resf,:) = nansem(resB{resf}.auc);
    for s = subjArray
        fn = sprintf('s%02d_auc',s);
        toR.(fn)(resf,:) = resB{resf}.auc(s);
    end
end

out=['name';toR.resname'];
fn = fieldnames(toR);
for f = 2:length(fn)
    vals = toR.(fn{f});
    assert(isfloat(vals));
    out = horzcat(out,[fn{f}; num2cell(vals)]);
end
 
if ~isempty(saveName)
    cell2csv(fullfile(resDir, ['aucSummary_' saveName '.csv']), out, ',', 2000);
end


