function [res, psyphys, idx, gRes] = CM_fMRIBehAnalysis(par,task)
% [res psyphys idx gRes] = CM_fMRIBehAnalysis(par)
% <par> = subject parameters setup by CM_Param

%specify variables for analysis
SUBJARRAY = [1,3:10,12:26];
EXCLUDEHIGHRTS = 0;
FIXDOUBLERESPS = 0;
%initialize variables
res = []; psyphys = []; idx = []; gRes = [];

if nargin < 2
    task = 'ret';
end

if ~isstruct(par)
    assert(ismember(par,SUBJARRAY));
    par = CM_Params(par, task);
end

%cd to subjects behavioral directory
%if ~exist(par.behavdir,'dir')
 %   mkdir(par.behavdir);
%end
%cd(par.behavdir);

%trialData compiles the data from the scans
dFN_test = dir(fullfile(par.logdir, '*test*.mat'));
fileNames_test = {dFN_test.name};
dotFiles_test = cell2mat(cellfun(@(x) x(1)=='.', fileNames_test, 'UniformOutput', false));
fileNames_test(find(dotFiles_test)) = []; %remove hidden files that are prepended with dots.
trialData = combineDataFile(fileNames_test, par.logdir);
numTrials = length(trialData);
trialRunNum = [];
for runNum = 1:par.(par.task).numRuns
    trialRunNum = [trialRunNum; zeros(par.(par.task).trialsPerRun(runNum),1)+runNum];
end

%reorder trialData so that explicit precedes countermeasures
trialDataOld = trialData;
for i =1:4
    trialData(i) = trialDataOld(4+i);
    trialData(4+i) = trialDataOld(i);
end

%fix some subject-specific idiosyncracies in logfiles
switch par.subNo
    case 1
        for i = 1:2:8
            trialData(i).SubHand='L';
            trialData(i+1).SubHand='R';
        end
    case 3
        %no index finger came through for run 1. grabbing from CM001_onsets_vector_test.xls
        trialData(1).resp = {'4$'	'3#'	'4$'	'4$'	'3#'	'4$'	'3#'	'3#'	'3#'	'4$'	'3#'	'3#'	'3#'	'4$'	'3#'	'3#'	'3#'	'4$'	'4$'	'4$'	'4$'	'4$'	'4$'	'3#'	'4$'	'3#'	'4$'	'3#'	'3#'	'3#'	'4$'	'4$'	'4$'	'3#'	'4$'	'3#'	'4$'	'4$'	'3#'	'3#'	'4$'	'4$'	'3#'	'4$'	'4$'	'3#'	'3#'	'3#'	'4$'	'3#'};     
        %something also wrong with run 2
        tmp = [8 8 9 9 9 8 9 8 8 9 8 8 9 8 8 8 9 8 8 8 8 9 9 8 8 8 8 9 8 9 9 9 9 9 9 8 8 8 9 9 8 9 8 8 8 8 8 8 8 9 ] ;
        trialData(2).resp = cell(size(tmp));
        trialData(2).resp(tmp==8) = {'8*'};
        trialData(2).resp(tmp==9) = {'9('};
    case 9
        %no data for run 1. grabbing from CM009_onsets_vector_test.xls
        trialData(1).resp = {'3#'	'4$'	'3#'	'3#'	'4$'	'4$'	'4$'	'3#'	'3#'	'4$'	'4$'	'4$'	'3#'	'4$'	'3#'	'4$'	'4$'	'3#'	'4$'	'4$'	'4$'	'3#'	'4$'	'4$'	'3#'	'4$'	'4$'	'4$'	'3#'	'3#'	'3#'	'4$'	'3#'	'3#'	'4$'	'4$'	'3#'	'4$'	'3#'	'4$'	'3#'	'4$'	'4$'	'3#'	'4$'	'4$'	'3#'	'3#'	'3#'	'4$'};
    case 17
        %for some reason run 8 has subtask == 1, which means nothing
        trialData(8).SubTask = 3;
end

for i = 1:8
    trialData(i).onset = (i-1)*(512)+4:10:i*(512)-18;
    trialData(i).hand = repmat(trialData(i).SubHand,1,50);
    trialData(i).task = repmat(trialData(i).SubTask,1,50);
end

switch par.task
    case 'ex'
        trialData = trialData(1:4);
    case 'cm'
        trialData = trialData(5:8);
end

resp = cat(2,trialData.resp);
idx.rt = cat(2,trialData.respRT)';
validTestTrials = ismember(resp,{'3#','4$','8*','9('})'; 

idx.excludedRts = [];
if EXCLUDEHIGHRTS
    idx.excludedRts = rt > par.rtThresh;
    validTestTrials = validTestTrials * ~excludedRts;
end

idx.time0 = cat(2,trialData.onset)';
idx.hand = cat(2,trialData.hand)';

idx.indexFinger = ismember(resp,{'4$', '9('})';
idx.middleFinger = ismember(resp,{'3#','8*'})';

idx.leftHand = ismember(resp,{'4$','3#'})';
idx.rightHand = ismember(resp,{'9(','8*'})';

idx.wrongHand = (idx.leftHand .* idx.hand == 'R') + (idx.rightHand .* idx.hand == 'L');
validTestTrials = validTestTrials .* ~idx.wrongHand;

idx.leftIndex = idx.indexFinger .* idx.leftHand;
idx.leftMiddle = idx.middleFinger .* idx.leftHand;
idx.rightIndex = idx.indexFinger .* idx.rightHand;
idx.rightMiddle = idx.middleFinger .* idx.rightHand;

oldNew = cat(1,trialData.cond);
idx.old = strcmp(oldNew , 'OLD');
idx.new = strcmp(oldNew , 'NEW');

idx.aa = strcmp(cat(1,trialData.group),'AA');
idx.ea = strcmp(cat(1,trialData.group),'EA');

idx.task = cat(2,trialData.task)';
idx.ex = (idx.task == 2);
idx.cm = (idx.task == 3);
idx.s_old = idx.indexFinger .* idx.ex + idx.middleFinger .* idx.cm;
idx.s_new = idx.indexFinger .* idx.cm + idx.middleFinger .* idx.ex;


idx.hit = idx.s_old .* idx.old;
idx.fa = idx.s_old .* idx.new;
idx.miss = idx.s_new .* idx.old;
idx.cr = idx.s_new .* idx.new;

idx.aaHit = idx.hit .* idx.aa;
idx.eaHit = idx.hit .* idx.ea;
idx.aaCr = idx.cr .* idx.aa;
idx.eaCr = idx.cr .* idx.ea;
idx.aaFa = idx.fa .* idx.aa;
idx.eaFa = idx.fa .* idx.ea;
idx.aaMiss = idx.miss .* idx.aa;
idx.eaMiss = idx.miss .* idx.ea;

idx.exLeftHand = idx.leftHand .* idx.ex;
idx.cmLeftHand = idx.leftHand .* idx.cm;
idx.exRightHand = idx.rightHand .* idx.ex;
idx.cmRightHand = idx.rightHand .* idx.cm;

idx.exLeftIndex = idx.exLeftHand .* idx.indexFinger;
idx.exLeftMiddle = idx.exLeftHand .* idx.middleFinger;
idx.exRightIndex = idx.exRightHand .* idx.indexFinger;
idx.exRightMiddle = idx.exRightHand .* idx.middleFinger;
idx.cmLeftIndex = idx.cmLeftHand .* idx.indexFinger;
idx.cmLeftMiddle = idx.cmLeftHand .* idx.middleFinger;
idx.cmRightIndex = idx.cmRightHand .* idx.indexFinger;
idx.cmRightMiddle = idx.cmRightHand .* idx.middleFinger;

idx.exHit = idx.ex .* idx.hit;
idx.cmHit = idx.cm .* idx.hit;
idx.exCr = idx.ex .* idx.cr;
idx.cmCr = idx.cm .* idx.cr;
idx.exFa = idx.ex .* idx.fa;
idx.cmFa = idx.cm .* idx.fa;
idx.exMiss = idx.ex .* idx.miss;
idx.cmMiss = idx.cm .* idx.miss;

idx.junk = ~(validTestTrials);

idx.allTrialOns = idx.time0;
% idx.allTrialOns = idx.time0 + ((idx.runNum-1) .* par.cam.numvols(idx.runNum)'* par.TR);

%% results sections
% res.pRespOld = sum(idx.respRFK)/sum(validTestTrials);
 res.pHit = sum(idx.hit)/sum(idx.old);
 res.pFa = sum(idx.fa)/sum(idx.new);
 res.pExHit = sum(idx.hit)/sum(idx.old .* idx.ex);
 res.pCmHit = sum(idx.hit)/sum(idx.old .* idx.cm);
 res.pExFa = sum(idx.fa)/sum(idx.old .* idx.ex);
 res.pCmFa = sum(idx.fa)/sum(idx.old .* idx.cm);
% 
% %dprime
 res.dprime = norminv(res.pHit) - norminv(res.pFa);
 res.exDprime = norminv(res.pExHit) - norminv(res.pExFa);
 res.cmDprime = norminv(res.pCmHit) - norminv(res.pCmFa);
% res.rtHit = nanmean(rt(find(idx.hit)),1);
% 
% %%psychophysics

% i=0;
% for phaseNum = 1:3
%     i = i+1;
%     thisPhase = (phase==phaseNum);
%     thisPhaseOld = thisPhase .* idx.old;
%     %1 hits by phase
%     psyphys.hitByPhz.mean(i) = sum(idx.hit .* thisPhaseOld)/sum(thisPhaseOld);
%     psyphys.hitByPhz.SE(i) = sqrt((psyphys.hitByPhz.mean(i) * (1-(psyphys.hitByPhz.mean(i))))/sum(thisPhaseOld));
%     psyphys.hitByPhz.ticks.x = 1:3;
%     psyphys.hitByPhz.ticks.y = [0 1];
%     psyphys.hitByPhz.xlabel = 'Phase ( 1 2 3)';
%     psyphys.hitByPhz.ylabel = 'Percent Hits';
%     psyphys.hitByPhz.xVals = 1:3;
%     psyphys.hitByPhz.marker = 'o-';
% 
% end


end
    
    



