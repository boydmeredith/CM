mask={ 'SEPT09_MVPA_MASK_resliced4mm.nii' 'rmotorMask_locerebellum_spmMaskminusSEPT09.img' 'rleftMtlMask_tbm112514.nii' 'rleftSplMask_tbm112514.nii' 'rleftAngMask_tbm112514.nii' 'rrightMtlMask_tbm112514.nii' 'rrightSplMask_tbm112514.nii'  'rrightAngMask_tbm112514.nii' 'rleftIfgAnat.nii', 'rleftBa44_45_47anat.nii'};
%classify late trs for all rois

for j = 1:length(mask)
    for scram = 0:1
        %late trs
        CM_run_mvpa_v4_tbm([],'','which_traintest',[1 2 3 4 1 2 3 4],'condNames',{'cmHit' 'cmCr'},'scramble', scram,'roiName',mask{j})
        classStr = {sprintf('*cmHit*%s*mat',mask{j}(1:10))}
    end
    cd('~/cm/Data/Functional/mvpa_results/');
    CM_batchMvpaPostProc(classStr,sprintf('cmtocm_late_%s',masks{j}(1:10)),pwd,0, {[]})
    for i = [1:6]
        testTrs = zeros(1,6);
        testTrs(i) = 1
        CM_run_mvpa_v4_tbm([],'','which_traintest',[1 2 3 4 1 2 3 4],'condNames',{'cmHit' 'cmCr'},'scramble', scram,'trWeights_train',testTrs, 'trWeights_test',testTrs,'roiName',mask{j})
    end
    classStr = {sprintf('*cmHit*1_*1_*%s*mat',mask{j}(1:10))}
    CM_batchMvpaPostProc(classStr,sprintf('cmtocm_trbytr_%s',masks{j}(1:10)),pwd,0, {[]})
    
    
end
