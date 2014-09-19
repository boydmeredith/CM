function CM_MakeRegs_MVPA(par, saveit)
%LTCam_MakeRegs(par, saveit) creates onsets and regresors (boxcar and motion)
%for use in SPM's GLM.
%inputs:
%par = subject-specific parameters created by LTCam_Params
%saveit = 1 if we'd like to save mat files for onsets and boxcar + motion
%regressors(default is 1)
%

analysis = 'mvpa';

if nargin<3
    saveit=1;
end

if length(par) > 1
    for i =1:length(par)
        CM_MakeRegs_MVPA(par(i))
    end
    return
end
        
if ~isstruct(par)
    par = CM_Params(par, 'ret',1);
end
if ~exist(par.behavdir,'dir')
    mkdir(par.behavdir);
end
cd (par.behavdir);

thisAnalysisDir = fullfile(par.subdir, analysis);
if ~exist(thisAnalysisDir,'dir')
    mkdir(thisAnalysisDir);
end

stimDur = 0;

switch par.task
    case 'ret'
        [~, ~, idx] = CM_fMRIBehAnalysis(par);
end

i = 0;

idxFnames = fieldnames(idx);
for xx  = 1:length(idxFnames)
    i = i+1;
    fName = idxFnames{xx};
    idx.thisOns = idx.(fName);
    stimOnsets{i} = idx.allTrialOns(find(idx.thisOns));
    stimNames{i} = fName;
    stimDurations{i} = stimDur;
end

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;

if ~exist(thisAnalysisDir)
    mkdir (thisAnalysisDir);
end

sessReg = zeros(sum(par.(par.task).numvols),length(par.(par.task).numvols) -1);
runStart = 0;
runEnd = 0;
for i =1:(length(par.(par.task).numvols) -1)
    runStart = runEnd + 1;
    runEnd = runStart+par.(par.task).numvols(i)-1;
    sessReg(runStart:runEnd,i) = ones(par.(par.task).numvols(i),1);
end

motRegs = [];
for runNum =1:8
    if runNum<5
        fid = fopen(fullfile(par.subdir, sprintf('TestExp/TestExp%02d/rp_aCM%03d_test%02d_0004.txt',runNum,par.subNo,runNum)));
    else
        fid = fopen(fullfile(par.subdir, sprintf('TestCM/TestCM%02d/rp_aCM%03d_CM%02d_0004.txt',runNum-4,par.subNo,runNum-4)));
    end
    motRegs = vertcat(motRegs,cell2mat(textscan(fid,'%n%n%n%n%n%n')));
end
fclose('all');

R = horzcat(sessReg,motRegs);

cd(thisAnalysisDir);
if saveit
    save mvpa_ons.mat onsets durations names;
    
%     save regs.mat R
end
