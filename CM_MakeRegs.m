function LTCam_MakeRegs(par, analysis saveit)
%LTCam_MakeRegs(par, saveit) creates onsets and regresors (boxcar and motion)
%for use in SPM's GLM. 
%inputs: 
%par = subject-specific parameters created by LTCam_Params
%saveit = 1 if we'd like to save mat files for onsets and boxcar + motion
%regressors(default is 1)
%



if nargin<2
    saveit=1;
end

if ~isstruct(par)
    par = LT_Params(par, 'cam',1);
end

cd (par.behavdir);

thisAnalysisDir = fullfile(par.analysisdir, analysis);

stimDur = 0;

switch par.task 
    case 'lab'
        [~, ~, idx] = LTLab_fMRIBehAnalysis_Ret(par);
    case 'cam'
        [~, ~, idx] = LTCam_fMRIBehAnalysis_Ret(par);
end

perfSet = {'hitSrcCor', 'hitSrcInc' 'hitNoSrc', 'miss', 'CR', 'srcFA', 'noSrcFA'};
srcSet = {'obj', 'face', 'scene', 'new'};
respSet = {'respObj', 'respFace', 'respScene', 'respOld' ,'respNew'};

i = 0;
for p = 1:length(perfSet)
    for srcNum = 1:length(srcSet)
        for respNum = 1:length(respSet)
            idx.thisOns = idx.(perfSet{p}) .* idx.(srcSet{srcNum}) .*  idx.(respSet{respNum});
            
            if ~isempty(idx.allTrialOns(find(idx.thisOns)))
                i = i+1;
                fName = sprintf('%s_%s_%s', perfSet{p}, srcSet{srcNum}, respSet{respNum});
                stimOnsets{i} = idx.allTrialOns(find(idx.thisOns));
                stimNames{i} = fName;
                stimDurations{i} = stimDur;
            end
        end
    end
end
              

if sum(idx.junk>0)
    i= i+1;
    stimOnsets{i} = idx.allTrialOns(find(idx.junk));
    stimNames{i} = 'junk';
    stimDurations{i} = 0;
end

onsets = stimOnsets;
names = stimNames;
durations = stimDurations;

if ~exist(thisAnalysisDir)
    mkdir (thisAnalysisDir);
end

sessReg = zeros(sum(par.lab.numvols),length(par.lab.numvols) -1);
for i =1:(length(par.lab.numvols) -1)
    sessReg((i-1)*par.lab.numvols(i)+1:i*par.lab.numvols(i),i) = ones(par.lab.numvols(i),1);
end

motRegs = [];
for runNum =1:length(par.(par.task).numvols)
    fid = fopen(fullfile(par.funcdir, sprintf('test%02d/rp_avol0005.txt',runNum)));
    motRegs = vertcat(motRegs,cell2mat(textscan(fid,'%n%n%n%n%n%n')));
end


R = horzcat(sessReg,motRegs);

cd(thisAnalysisDir);
if saveit
    save ons.mat onsets durations names;    
    save regs.mat R
end