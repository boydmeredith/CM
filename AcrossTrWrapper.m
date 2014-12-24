mask = {'rmotorMask_locerebellum_spmMaskminusSEPT09.img' 'SEPT09_MVPA_MASK_resliced4mm.nii' 'rleftMtlMask_tbm112514.nii' 'rleftSplMask_tbm112514.nii' 'rleftAngMask_tbm112514.nii' 'rrightMtlMask_tbm112514.nii' 'rrightSplMask_tbm112514.nii'  'rrightAngMask_tbm112514.nii'};
%should probably repeat with 'rleftIfgAnat.nii', 'rleftBa44_45_47anat.nii'
masksToClassify = mask(6:end)
%classify late trs for all remaining rois
for scram = 0:1
   for j = 1:length(masksToClassify)
        
        %all rois for late trs
        CM_run_mvpa_v4_tbm([],'','which_traintest',[1 2 3 4 1 2 3 4 ],'condNames',{'exHit' 'exCr'},'scramble', scram,'roiName',masksToClassify{j})
        CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 1 1 2 2 2 2],'condNames',{'hit' 'cr'},'scramble',scram,'roiName',masksToClassify{j})
        
    end
end
cd('~/cm/Data/Functional/mvpa_results/')
%all late rois
CM_batchMvpaPostProc({'*ex*ex*1_2_3_4_1_2_3_4*33*mat' '*hitVcr*1_1_1_1_2_2_2_2*33*mat'}, 'allrois_latetrs',pwd,0, {[], 2})

mask = {'rmotorMask_locerebellum_spmMaskminusSEPT09.img' 'SEPT09_MVPA_MASK_resliced4mm.nii' 'rleftMtlMask_tbm112514.nii' 'rleftSplMask_tbm112514.nii' 'rleftAngMask_tbm112514.nii' 'rrightMtlMask_tbm112514.nii' 'rrightSplMask_tbm112514.nii'  'rrightAngMask_tbm112514.nii'};

%classify trs 3:6 for whole brain roi and motor roi
masksToClassify= mask(1:2)

for scram = 0:1      
     for i = [3:6]
         for j = 1:length(masksToClassify)
             testTrs = zeros(1,6);
             testTrs(i) = 1
             CM_run_mvpa_v4_tbm([],'','which_traintest',[1 2 3 4 1 2 3 4 ],'condNames',{'exHit' 'exCr'},'trWeights_train',testTrs, 'trWeights_test',testTrs,'scramble', scram,'roiName',masksToClassify{j})
             if i > 3
                CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 1 1 2 2 2 2],'condNames',{'hit' 'cr'},'trWeights_train',testTrs, 'trWeights_test',testTrs,'scramble',scram,'roiName',masksToClassify{j})
             end
         end
     end
end

%process it all
cd('~/cm/Data/Functional/mvpa_results/')
%trbytr whole brain and motor
CM_batchMvpaPostProc({'*ex*ex*1_2_3_4_1_2_3_4*locer*mat' '*ex*ex*1_2_3_4_1_2_3_4*SEPT*mat' '*hitVcr*1_1_1_1_2_2_2_2*locer*mat' '*hitVcr*1_1_1_1_2_2_2_2*SEPT*mat'}, 'wb_motor_trbytr',pwd,0, {[],[],2, 2})
