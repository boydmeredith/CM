function [resB cvB out] = CM_batchMvpaPostProc(classStr, saveName, resDir, plotit, xvalIterToReport, runs_to_report,subjArray, task)
%function [resB cvB out] = CM_batchMvpaPostProc(classStr, saveName, resDir, plotit, xvalIterToReport, runs_to_report,subjArray, task)
%inputs:
%<resDir> the directory containing the results mat files
%<classSt>: default value is '*conds*". Filters mat files
%<xvalIterToReport>: allows specification of an array of xval iterations to report. This is specifically for cases when the xval iterations are very different, such as when one xval bin is all explicit trials and one xval bin is all countermeasures trials
%<task>: default value is 'ret'

toR = [];

if isempty(classStr)
    classStr = '*conds*';
end
if isempty(saveName)
    saveName = classStr;
end

if nargin < 3 
    resDir = pwd
end
if nargin < 4
    plotit = 0
end
if nargin < 5
    xvalIterToReport = []
end
if nargin < 6
    runs_to_report = []
end
if nargin < 7
    subjArray = []
end
if nargin < 8
    task = 'ret';
end

userSubjArray = subjArray;

% get list of all files matching classStr and make sure there's something
% there
fileNames = '';
for iclassStr = 1:length(classStr)
	dFN = dir(fullfile(resDir, [classStr{iclassStr}]));
	fileNames = horzcat(fileNames, {dFN.name});
end
dotFiles = cell2mat(cellfun(@(x) x(1)=='.', fileNames, 'UniformOutput', false));
fileNames(find(dotFiles)) = []; %remove hidden files that are prepended with dots.
nfile = length(fileNames)
assert(nfile>0);

% iterate through the files and get results
for ifile =1:nfile
    % tell the user what we're working on
    display(['Now processing ' fileNames{ifile} '...']);
    % load results
    thisRes =load(fileNames{ifile},'-mat');
    % If it is an across subject classification, get the res structure into
    % the epected format
    isAcrossSubs = isfield(thisRes,'results');
    res.isAcrossSubsResStruct = isAcrossSubs;
    if isAcrossSubs
        res.subj{1}.penalty.nVox.weights.iter{1} = thisRes.results{1};
        res.subj{1}.penalty.nVox.weights.expt{1}.saveName = fileNames{ifile};
        res.subj{1}.penalty.nVox.weights.condensed_runs = thisRes.results{1}.condensed_runs;
        subjArray = 1;
    else
        res = thisRes.res;
    end
    
    % if user did not specify subject array, use info in res struct
    if isempty(userSubjArray)
        subjArray = res.subjArray;
    end
    % do the actual processing of this results struct
    [resB{ifile} cvB{ifile}] = CM_mvpaPostProc(res, xvalIterToReport,runs_to_report, task, plotit, fullfile(resDir,'figs'), subjArray);
    
    % package the processed results with the name, subject array, and expt
    % params used
    resB{ifile}.name = fileNames{ifile};
    cvB{ifile}.name = fileNames{ifile};
    resB{ifile}.subjArray = subjArray;
    resB{ifile}.expt = res.subj{1}.penalty.nVox.weights.expt{1};

end

% now that we have a structure of results, we can put iterate through them
% and input them into a data frame format to send to R
fieldsToRecord = {'condNames','which_traintest','roiName','num_results_iter','scramble','trWeights','trWeights_train','trWeights_test'};
nres = length(resB); % number of results to record
idxToEdit = 0; % counter for index in data frame

%iterate through each result
for ires = 1:nres
	subArray = resB{ires}.subjArray;
	nsubs = length(resB{ires}.subjArray);
    %iterate through each subject 
	for jsub = subArray
		idxToEdit = idxToEdit+1;
		s_ix = find(ismember(subArray,jsub));
        % record result values of interest
		toR.name{idxToEdit} = resB{ires}.name;
		toR.mean_auc(idxToEdit) = nanmean(resB{ires}.auc);
		toR.sem_auc(idxToEdit) = nansem(resB{ires}.auc);
		toR.subNo{idxToEdit} = sprintf('s%02d',jsub);
		toR.auc(idxToEdit) = resB{ires}.auc(jsub);
        toR.xvalIterReported{idxToEdit} = num2str(xvalIterToReport);
        
        %record values of interest from expt params
        for kfield = 1:length(fieldsToRecord)
            fieldname = fieldsToRecord{kfield}
            if isfield(resB{ires}.expt, fieldname)
                %use strjoin to turn value into a delimited string
                value = strjoin(resB{ires}.expt.(fieldname),';');
            else
                value='';
            end
            toR.(fieldname){idxToEdit} = value;
        end
    end
end
 
% take the toR struct and turn into a big matrix that can actually be
% turned into a csv
 
fn = fieldnames(toR);
out=['name';toR.name'];
fn(strcmp(fn,'name'))=[];

for f = 1:length(fn)
    vals = toR.(fn{f});
    if(isfloat(vals))
	    out = horzcat(out,[fn{f}; num2cell(vals)']);
    else
	    out = horzcat(out,[fn{f}; vals']);
    end

end

%turn the big matrix into a csv

if ~isempty(saveName)
    cell2csv(fullfile(resDir, ['aucSummary_' runs_to_report '_' saveName '.csv']), out, ',', 2000);
end

end
