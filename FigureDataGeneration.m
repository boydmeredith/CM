NSUBS = 24;
subj_array  = [1,3:10,12:26];
%generate data for across subject hit v cr classification masked by ex hit
%> cr and cm hit > cr effects inclusively masked at the 2nd level across
%the held-out 24 subjects.
TRSTOTEST = 1:6;
for i =TRSTOTEST
    display(['Now running across subject classification for TR ' num2str(i) '... ']);
    trW = [0 0 0 0 0 0];
    trW(i) = 1;
    CM_wrap_mvpa_across_subjects_masks('EXminus',[],[], 'trWeights',trW);
end
%repeat for late TRs
CM_wrap_mvpa_across_subjects_masks('EXminus',[],[]);

%results will likely require that I sort out where the condensed runs live
%in the resulting res results struct
cd('/Users/Jesse/fMRI/COUNTERMEASURES/Data/Functional/mvpa_results');
%examine results for all 144
theseResStr = '*trEXteAll*trW*minus*'
assert(length(dir(theseResStr)) == (1+length(TRSTOTEST))*NSUBS);
%analyze for CM runs
CM_batchMvpaPostProc(pwd,theseResStr,2,'ret',[],'acrossSs_trEXteCM_combinedFXMask_hitvcr',0,[5:8]);
%analyze for EX runs
CM_batchMvpaPostProc(pwd,theseResStr,2,'ret',[],'acrossSs_trEXteEX_combinedFXMask_hitvcr',0,[1:4]);

% generate
for iTrainTr = 3:4
    for jTestTr = 1:6
        trW_train = zeros(1,6);
        trW_test = zeros(1,6);
        trW_train(iTrainTr) = 1;
        trW_test(jTestTr) = 1;
        
        CM_run_and_report_mvpa(subj_array,[],2,'trWeights_train',trW_train,...
            'trWeights_test',trW_test, 'which_traintest',[1 1 1 1 2 2 2 2]);
        
        CM_run_and_report_mvpa(subj_array,[],[],'trWeights_train',trW_train,...
            'trWeights_test',trW_test, 'which_traintest',[1 1 2 2 0 0 0 0]);

    end
end