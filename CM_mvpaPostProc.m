function [res, cv, idx] = CM_mvpaPostProc(qqq, xvalIterToReport,runs_to_report,task, plotit, saveDir, subjArray)
%function [res, cv, idx] = CM_mvpaPostProc(qqq, xvalIterToReport, runs_to_report,task, plotit, saveDir, subjArray)

auc = [];
if isempty(subjArray)
    subjArray = qqq.subjArray;
end

saveTag = sprintf('xvals%s_runsrpotred%s_',strrep(num2str(xvalIterToReport,' ','')),strrep(num2str(runs_to_report,' ','')));

qqq = zeroUnwantedTrials(qqq,subjArray,xvalIterToReport,runs_to_report)
for s = subjArray
    [cv resS] = CM_singleSubjProc(s, qqq, xvalIterToReport, task);
    
    fn = fieldnames(cv.reordered);
    for f = 1:length(fn)
        toR.(fn{f}){s} = asColumn(double(cv.reordered.(fn{f})));
    end
    
    toR.subs{s} = s*ones(size(toR.(fn{f}){s}));
end
%% write out data across subjects
fn = fieldnames(toR);

for f = 1:length(fn)
    out.(fn{f}) = vertcat(toR.(fn{f}){:});
    toPrint(:,f) = [fn{f}; num2cell(out.(fn{f}))];
end

if ~isempty(saveDir)
	cell2csv(fullfile(saveDir, [saveTag resS{1}.saveName '.csv']), toPrint, ',', 2000);
end

res.auc = nan(size(subjArray));
res.far =  nan(length(subjArray),80);
res.hr =  nan(length(subjArray),80);
for s=subjArray
	ix = out.subs==s;
	[res.auc(s), res.far(s,:), res.hr(s,:)]= getAuc(out,ix,0);
end

if plotit
   plot_roc_logits(res,out)
   if ~isempty(saveDir)
        figurewrite([saveTag resS{1}.saveName],[],[],[saveDir],[]);
    end
end
end

function plot_roc_logits(res, out)
    subplot(1,3,1);
    plot(nanmean(res.far),nanmean(res.hr));
    errorbar3(nanmean(res.far),nanmean(res.hr),ste(res.hr),1,[.2 .6 1]);
    title(sprintf('AUC = %.2f',nanmean(res.auc)));
    axis square;
    subplot(1,3,2);
    hist(out.actsVec( out.desiredsVec==1),25);
    title('class A probs for A trials');
    axis square;
    subplot(1,3,3);
    hist(out.actsVec( out.desiredsVec==2),25);
    title('class A probs for B trials');
    axis square;
end


function [cv resS] = CM_singleSubjProc(s, qqq, xvalIterToReport, task)
% behavioral analysis
par = CM_Params(s, task, 0);
[~,~,idxB] = CM_fMRIBehAnalysis(par);
resS = qqq.subj{s}.penalty.nVox.weights.expt;
for resit = 1:length(qqq.subj{s}.penalty.nVox.weights.iter)
    % read in classification data
    s, resit
    if xvalIterToReport
        resStruct = qqq.subj{s}.penalty.nVox.weights.iter{resit}.iterations(xvalIterToReport);
    else
        resStruct = qqq.subj{s}.penalty.nVox.weights.iter{resit}.iterations;
    end
    
    %allOns = sort([resS.onsets_test_in_classifier{:}]);
%     testOns = zeros(size(idxB.allTrialOns));
%     for i = 1:length(resS{1}.condNames)
%         testOns = testOns+idxB.(resS{1}.condNames{i});
%     end
%     inClassifier = find(testOns);
    
    %global signal stuff
    %globalSignalDat = load('/biac4/wagner/biac3/wagner5/alan/perceptMnemonic/fmri_data/mvpa_files/pre2013/MeanSigIntensityInOTCortex');
    %globalSig = globalSignalDat.res.subj{s}.penalty.nVox.weights.iter{1}.iterations.acts(1,:);
    
    %% extract raw data from the resStruct
    for i = 1:length(resStruct)
        
        cv.raw.correctsVec{resit*i} = resStruct(i).perfmet.corrects; %#ok<*AGROW> %correct
        cv.raw.actsVec{resit*i} = resStruct(i).acts(1,:) ; % evidence for class 1
        cv.raw.actsVec2{resit*i} = resStruct(i).acts(2,:) ; % evidence for class2
        
        cv.raw.desiredsVec{resit*i} = resStruct(i).perfmet.desireds; % correct class
        
        cv.raw.actsOfTrueClass{resit*i} = cv.raw.actsVec{resit*i}; % evidence for true class
        cv.raw.actsOfTrueClass{resit*i}(cv.raw.desiredsVec{resit*i}~=1) = cv.raw.actsVec2{resit*i}(cv.raw.desiredsVec{resit*i}~=1);
        
        cv.raw.actsOfFalseClass{resit*i} = cv.raw.actsVec2{resit*i}; % evidence for false class
        cv.raw.actsOfFalseClass{resit*i}(cv.raw.desiredsVec{resit*i}~=1) = cv.raw.actsVec{resit*i}(cv.raw.desiredsVec{resit*i}~=1);
        
        cv.raw.actsDifference{resit*i} = cv.raw.actsVec{resit*i} - cv.raw.actsVec2{resit*i}; % difference in evidence for true vs. false class
        cv.raw.actsDifference{resit*i}(cv.raw.desiredsVec{resit*i}~=1) = -1*cv.raw.actsDifference{resit*i}(cv.raw.desiredsVec{resit*i}~=1);
        
        cv.raw.guessesVec{resit*i} = resStruct(i).perfmet.guesses; % binary classifier guesses
        cv.raw.signedActs{resit*i} = abs(cv.raw.actsVec{resit*i} - .5) .* (2*cv.raw.correctsVec{resit*i} - 1); % acts signed in correct direction.
        
        cv.raw.xval_iter{resit*i} = i+zeros(size(cv.raw.correctsVec{resit*i}));
        
        if isfield(resStruct(i), 'test_idx')
            testIdx{resit*i} = resStruct(i).test_idx;
        else
            testIdx{resit*i} = 1:length(cv.raw.guessesVec{resit*i});
        end
        
        cv.raw.signedEv{resit*i} = logit(cv.raw.actsVec{resit*i}); %take the logit of the probabilistic estimate for class 1
        cv.raw.unsignedEv{resit*i} = logit(cv.raw.actsOfTrueClass{resit*i});%take the logit of the probabilistic estimate for the correct class
        
        nTrainingTrials(resit*i) = length(resStruct(i).train_idx);
    end
    
    testIdxCat = [testIdx{:}];
    [~, testIdxInv] = sort(testIdxCat);
    
    fn = fieldnames(cv.raw);
    
    % consoildate vectors across classifier iterations
    % reorder them to match chronological presentation time.
    for f = 1:length(fn)
        cv.consolidated.(fn{f}) = [cv.raw.(fn{f}){:}];
        cv.reordered.(fn{f}) = cv.consolidated.(fn{f})(testIdxInv);
    end
    
    %compare onsets from qqq to onsets from idxB
    %onsInClassifier = ismember(idxB.alltrials, allOns);
    
    %% get behavioral data
    %     if strcmp(task, 'mnem')
    %         idx.cor = idxB.cor(onsInClassifier);
    %     elseif strcmp(task, 'perc')
    %         idx.cor = idxB.cor(onsInClassifier);
    %     end
    %
    %      idx.hit = idxB.hit(testIdxCat);
    
    
    %     dat{s}.idx = idx;
    %     dat{s}.cv = cv;
    %     datB{s} = idxB;
    
    %% results summary
    %     assert(nanmean(cv.reordered.desiredsVec(idx.(resS{1}.condNames{1}))) == 1);
    %     assert(nanmean(cv.reordered.desiredsVec(idx.(resS{1}.condNames{2}))) == 2);
    res.cor.class1(s) = nanmean(cv.reordered.correctsVec(cv.reordered.desiredsVec==1));
    res.cor.class2(s) = nanmean(cv.reordered.correctsVec(cv.reordered.desiredsVec==2));
    %     res.cor.class1Cor(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.cor'==1) ));
    %     res.cor.class2Cor(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.cor'==1) ));
    %     res.cor.class1Inc(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==1) & (idx.inc'==1) ));
    %     res.cor.class2Inc(s) = nanmean(cv.reordered.correctsVec((cv.reordered.desiredsVec==2) & (idx.inc'==1) ));
    %
    res.MeanNTrainingTrials(s) = nanmean(nTrainingTrials);
    res.sub(s).S = resS;

    end
end


function [auc_80_bins, fas_80, hits_80] = getAuc(out,ix,plotit)
% sort by signed classifier "confidence" (for ROI curves)
actsVec = out.actsVec(ix);
correctsVec = out.correctsVec(ix);
desiredsVec = out.desiredsVec(ix);
[sorted_diffs ind] = sort(actsVec,'descend');
correct_sorted = correctsVec(ind);
desireds_sorted = desiredsVec(ind);

for i = 1:length(sorted_diffs);
    hit_rate(i) = length(intersect(find(desireds_sorted == 1),[1:i])) / length(find(desireds_sorted == 1));
    fa_rate(i) = length(intersect(find(desireds_sorted == 2),[1:i])) / length(find(desireds_sorted == 2));
end
auc = auroc(hit_rate',fa_rate');

% create ROC function with 80 discrete bins (this allows the ROC curves to be more easily averaged across individuals)
roc_bin_intervals = .975:-.025:-1;
for bin_num = 1:80
    hits_80(bin_num)=length(intersect(find(desireds_sorted == 1),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 1));
    fas_80(bin_num)=length(intersect(find(desireds_sorted == 2),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 2));
end
auc_80_bins = auroc(hits_80',fas_80');

if plotit    
    plot(fa_rate,hit_rate,'.-');
    hold on;
    plot([0 1],[0 1],'r');
    set(gcf,'Color','w');
    %     xlabel('P(Old|New)')
    %     ylabel('P(Old|Old)')
end
end

function [out] = logit(prob)

out = log(prob ./ (1 - prob));

end

function [out] = asColumn(vec)
sz=size(vec);
nonSingletonDims = sz>1;

if sum(nonSingletonDims>1)
    error('data must only have one non-singleton dimension')
end

vecLength = sz(nonSingletonDims);

out = reshape(vec,vecLength, 1);

end

function res = zeroUnwantedTrials(res,subjArray,xvalIterToReport,runs_to_report)
for s = subjArray
    thisS = res.subj{s}.penalty.nVox.weights.iter;
    for r_it = 1:length(thisS)
        for xval = xvalIterToReport
            testidx = thisS{r_it}.iterations(xval).test_idx;
            switch runs_to_report
                case 'reportEX'
                    runs_to_report = [1:4];
                case 'reportCM'
                    runs_to_report = [5:8];
            end
            trials_to_report = testidx(ismember(thisS{r_it}.condensed_runs(testidx),runs_to_report))
            %zero out values that we don't want
            idx_to_report = find(ismember(testidx,trials_to_report));
            perfmet = thisS{r_it}.iterations(xval).perfmet;
            res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xval).perfmet.guesses = perfmet.guesses(idx_to_report);
            res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xval).perfmet.desireds = perfmet.desireds(idx_to_report);
            res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xval).perfmet.corrects = perfmet.corrects(idx_to_report);
            res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xval).acts=res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xvalIterToReport).acts(:,idx_to_report);
            res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xval).test_idx=res.subj{s}.penalty.nVox.weights.iter{r_it}.iterations(xval).test_idx(idx_to_report);
        end
    end
end

