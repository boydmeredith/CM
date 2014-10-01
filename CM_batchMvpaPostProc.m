function [resB cvB out] = CM_batchMvpaPostProc(resDir, classStr, xvalIterToReport, task, subjArray, saveName, plotit, restrictAcrossSubToTheseRuns) 
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

dFN = dir(fullfile(resDir, [classStr]));

fileNames = {dFN.name};
dotFiles = cell2mat(cellfun(@(x) x(1)=='.', fileNames, 'UniformOutput', false));
fileNames(find(dotFiles)) = []; %remove hidden files that are prepended with dots.

for resf =1:length(fileNames)
    thisRes =load(fileNames{resf},'-mat');
%     if length(res.subjArray) > length(res.subj)
%         display(['skipping over ' res.subj{1}.penalty.nVox.weights.expt{1}.saveName ' b/c the results struct is shorter than the subject array...\n']);
%         continue
%     end

    isAcrossSubsResStruct = isfield(thisRes,'results');
    if isAcrossSubsResStruct
        res.subj{1}.penalty.nVox.weights.iter{1} = thisRes.results{1};
        res.subj{1}.penalty.nVox.weights.expt{1}.saveName = fileNames{resf}; 
        subjArray = 1;
        if ~isempty(restrictAcrossSubToTheseRuns)
            testidx = res.subj{1}.penalty.nVox.weights.iter{1}.iterations(2).test_idx;
            if strcmp(restrictAcrossSubToTheseRuns,'reportEX');
                trials_to_report = testidx(ismember(res.subj{1}.penalty.nVox.weights.iter{1}.condensed_runs(testidx),1:4))
            elseif strcmp(restrictAcrossSubToTheseRuns,'reportCM')
                trials_to_report = testidx(ismember(res.subj{1}.penalty.nVox.weights.iter{1}.condensed_runs(testidx),5:8))
            else
                error('You supplied a bad string for restricting the analysis');
            end
            %zero out values that we don't want
            idx_to_report = find(ismember(testidx,trials_to_report));
            for xval = xvalIterToReport
                perfmet = res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).perfmet;
                res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).perfmet.guesses = perfmet.guesses(idx_to_report);
                res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).perfmet.desireds = perfmet.desireds(idx_to_report);
                res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).perfmet.corrects = perfmet.corrects(idx_to_report);
                res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).acts=res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xvalIterToReport).acts(:,idx_to_report)
                res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).test_idx=res.subj{1}.penalty.nVox.weights.iter{1}.iterations(xval).test_idx(idx_to_report)
            end
        end
    else
        res = thisRes.res;
    end
    [resB{resf} cvB{resf}] = CM_mvpaPostProc(res, xvalIterToReport, task, plotit, fullfile(resDir,'figs'), subjArray);
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
out=horzcat(['xval iterations reported'; repmat({num2str(xvalIterToReport)},resf,1)], out);

 
if ~isempty(saveName)
    cell2csv(fullfile(resDir, ['aucSummary_' restrictAcrossSubToTheseRuns '_' saveName '.csv']), out, ',', 2000);
end


