function [resB cvB] = batchMvpaPostProc(resDir, classStr, xvalIterToReport, task, subjArray) 
%function [resB cvB] = batchMvpaPostProc(dir, classStr, xvalIterToReport, task) 
%inputs:
%<resDir> the directory containing the results mat files
%<classSt>: default value is '*conds*". Filters mat files
%<xvalIterToReport>: allows specification of an array of xval iterations to report. This is specifically for cases when the xval iterations are very different, such as when one xval bin is all explicit trials and one xval bin is all countermeasures trials
%<task>: default value is 'ret'


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

for i =1:length(fileNames)
    load(fileNames{i});
    if length(res.subjArray) > length(res.subj)
        display(['skipping over ' res.subj{1}.penalty.nVox.weights.expt{1}.saveName ' b/c the results struct is shorter than the subject array...\n']);
        continue
    end
    [resB{i} cvB{i}] = CM_mvpaPostProc(res, xvalIterToReport, task, fullfile(resDir,'figs'), subjArray);
    resB{i}.name = fileNames{i};
    cvB{i}.name = fileNames{i};
end
    
