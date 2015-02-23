function par = CM_Params(subNo, task, loadScans, proclus)
% sets up parameters for batching of CM experiment.
% <subNo>: a vlue in the range of 1,3:10,12:26, or the subId (e.g. 'CM001')
% <load scans>:
% 1 - load all the scans (this is necessary for analyses of fmri data).
% 0 - do not load scans (saves time)

% this script is adapted by TBM from AG's PM_Params

%% basic task information
par.subNo = subNo;
if (nargin >= 2)
    par.task = task;
else
    par.task = 'ret';
end

%% convert input
if (nargin > 3)
    par.loadScans = loadScans;
else
    par.loadScans = 1; % by deault, load scans
end
if nargin < 4
    par.proclus = 0;
end

par.str = sprintf('CM%03d',par.subNo);

%% subject-specific information

%par.flagIt = should we flag this subject for some special reason?

par.ret.numvols = [256 256 256 256 256 256 256 256];
par.ret.scan0 = [4 4 4 4 4 4 4 4];
par.ret.trialsPerRun = [50 50 50 50 50 50 50 50];
par.ret.numRuns = 8;
par.ex.numvols = [256 256 256 256];
par.ex.scan0 = [4 4 4 4];
par.ex.trialsPerRun = [50 50 50 50];
par.ex.numRuns = 4;
par.cm.numvols = [256 256 256 256];
par.cm.scan0 = [4 4 4 4];
par.cm.trialsPerRun = [50 50 50 50];
par.cm.numRuns = 4;
par.TR = 2;
race = {'EA','','EA', 'AA','EA','EA','EA','AA','EA','AA','','EA','AA','EA','EA','EA','EA','EA','AA','AA','AA','EA','EA','EA','AA','EA'};
par.race = race{par.subNo};
switch par.subNo
    case 1
        %didn't do instructed countermeasures. Instead, he thought 'I don't
        %know that face' for recognized faces
        par.flagIt = 1; %don't expect an angular effect here
    case 8
        %wasn't sure to do with faces he was unsure about
    case 12
        %sleepy at encoding
        %AT THIS POINT INSTRUX CHANGE: suggested to subjects that they
        %perform encoding strategies during ISI
    case 14
        %sleepy at encoding (reported on all runs)
    case 22
        %tapping feet and moving hands
    case 25
        %sleepy at encoding AND retrieval towards end of each run
    case 26 
        %subject kept talking between scans to ask clarifying Qs
        %didn't fixate during final session
    
end

par.substr = par.str;
subject = par.str;
par.(par.task).numscans = sum(par.(par.task).numvols);


par.scansSelect.(par.task) = 1:length(par.(par.task).trialsPerRun);

%% select which scans to include
% fill in at a later date

%% directory information
if proclus
	par.exptdir = '/hsgs/projects/awagner/jtboyd/CM_ret/';
else
	par.exptdir = '/Users/Jesse/fMRI/COUNTERMEASURES';
end
par.scriptsdir = fullfile(par.exptdir,'Scripts');
par.datadir = fullfile(par.exptdir,'Data');
par.funcdir = fullfile(par.datadir,'Functional');
par.subdir = fullfile(par.funcdir, subject);
par.logdir = fullfile(par.subdir,'logfiles');
par.analysisdir = fullfile(par.subdir, 'analysis');
par.behavdir = fullfile(par.subdir, 'behav');


%univariate SPM setup
par.timing.fmri_t = 36;
par.timing.fmri_t0 = 35;
par.timing.units = 'secs';
par.volt = 1;
par.bases.hrf.derivs = [1 1];
par.cvi = 'AR(1)'; %note that this actually gets changed to AR(0.2) in spm_fmri_spm_ui.  
par.global = 'None';
par.mask = '';
par.constat = 'T';
par.sess.condstr = 'ons.mat';
par.sess.regstr = 'regs.mat';
par.sess.hpf = 128;

warning('off');
if ~isempty(par.substr) && par.loadScans
    [par.swascanfiles.(par.task) par.swascanfilesByRun.(par.task)] = findScans(par, 'swaCM*.nii',3,par.subdir);
else
end
warning('on');
end

function [scanfiles scansByRun] = findScans(par, matchstring, scanDim, scandir)
% input:
% <par> - the relevant par structure
% <matchstring> - a string representing a class of image data of interest, e.g. 'Rrun*.nii'
% <scanDim> - how many dimensions is the data? (e.g. 3)
% <scanDir> - where the functional data are located

% returns: 
% <scanfiles> - lists all the images of interest, for the specified
% <scansByRun> - lists all scan files of interest by run
scanfiles = [];
scansByRun = [];

clear sf_h sf_h2

if scanDim == 3
    for rnum=1:length(par.(par.task).numvols)
        if strcmp(par.task,'ex') || rnum<5
            runDir = fullfile(scandir, sprintf('TestExp/TestExp%02d',rnum));
        elseif strcmp(par.task,'cm') || rnum>=5
            runDir = fullfile(scandir, sprintf('TestCM/TestCM%02d',rnum-4));  
        end
        dFN = dir(fullfile(runDir, matchstring));
        scansByRun{rnum} = {dFN.name}';
        scansByRun{rnum} = strcat([runDir '/'],scansByRun{rnum});
        scanfiles = cat(1,scanfiles, scansByRun{rnum});
        
    end
elseif scanDim == 4
    error('we can''t handle 4d nifitis yet');
end
end

    
