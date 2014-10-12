NSUBS = 24;
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
