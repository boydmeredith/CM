function [resB cvB out] = CM_batchMvpaPostProc(resDir, classStr, xvalIterToReport, task, subjArray, saveName, plotit, runs_to_report)
%function [resB cvB out] = CM_batchMvpaPostProc(resDir, classStr, xvalIterToReport, task, subjArray, saveName, plotit, runs_to_report)
%inputs:
%<resDir> the directory containing the results mat files
%<classSt>: default value is '*conds*". Filters mat files
%<xvalIterToReport>: allows specification of an array of xval iterations to report. This is specifically for cases when the xval iterations are very different, such as when one xval bin is all explicit trials and one xval bin is all countermeasures trials
%<task>: default value is 'ret'

toR = [];
allSubsArray = [1 3:10 12:26];
userSubjArray = subjArray;


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

for ires =1:length(fileNames)
    display(['Now processing ' fileNames{ires} '...']);
    thisRes =load(fileNames{ires},'-mat');
    isAcrossSubsResStruct = isfield(thisRes,'results');
    if isAcrossSubsResStruct
        res.subj{1}.penalty.nVox.weights.iter{1} = thisRes.results{1};
        res.subj{1}.penalty.nVox.weights.expt{1}.saveName = fileNames{ires};
        res.subj{1}.penalty.nVox.weights.condensed_runs = thisRes.results{1}.condensed_runs;
        subjArray = 1;
    else
        res = thisRes.res;
    end
    res.isAcrossSubsResStruct = isAcrossSubsResStruct;
    if isempty(userSubjArray)
        subjArray = res.subjArray;
    end
    [resB{ires} cvB{ires}] = CM_mvpaPostProc(res, xvalIterToReport,runs_to_report, task, plotit, fullfile(resDir,'figs'), subjArray);
    resB{ires}.name = fileNames{ires};
    cvB{ires}.name = fileNames{ires};
    resB{ires}.subjArray = subjArray;
end

nres = length(resB);
idxToEdit = 0;
for ires = 1:nres
	nsubs = length(resB{ires}.subjArray);
	for jsub = 1:nsubs
		idxToEdit = idxToEdit+1;
		toR.resname{idxToEdit} = resB{ires}.name;
		toR.mean_auc(idxToEdit) = nanmean(resB{ires}.auc);
		toR.sem_auc(idxToEdit) = nansem(resB{ires}.auc);
		toR.subNo{idxToEdit} = sprintf('s%02',jsub);
		if jsub > length(resB{ires}.auc)
			toR.auc(idxToEdit);
		end
		toR.auc(ires,:) = resB{ires}.auc(jsub);
         end
 end
 
fn = fieldnames(toR);

out=['name';toR.resname'];

for f = 2:length(fn)
    vals = toR.(fn{f});
    if(isfloat(vals))
	    out = horzcat(out,[fn{f}; num2cell(vals)']);
    else
	    out = horzcat(out,[fn{f}; vals']);
    end

end
out=horzcat(['xval iterations reported'; repmat({num2str(xvalIterToReport)},idxToEdit,1)], out);

if ~isempty(saveName)
    cell2csv(fullfile(resDir, ['aucSummary_' runs_to_report '_' saveName '.csv']), out, ',', 2000);
end

end
