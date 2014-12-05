mask = {'rleftMtlMask_tbm112514.nii' 'rleftSplMask_tbm112514.nii' 'rleftAngMask_tbm112514.nii' 'rightMtlMask_tbm112514.nii' 'rrightSplMask_tbm112514.nii'  'rrightAngMask_tbm112514.nii'};
mask = { 'rrightMtlMask_tbm112514.nii' 'rrightSplMask_tbm112514.nii'  'rrightAngMask_tbm112514.nii'};

for j = 1:length(mask)
     for i = [1:6]
         for scram = 0:1
             testTrs = zeros(1,6);
             testTrs(i) = 1
             %CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 2 2 1 1 2 2],'condNames',{'exHit' 'exCr'},'trWeights_train',[0 0 1 0 0 0], 'trWeights_test',testTrs,'num_results_iter',1)
             CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 2 2 1 1 2 2],'condNames',{'exHit' 'exCr'},'trWeights_train',[0 0 1 0 0 0], 'trWeights_test',testTrs,'roiName',mask{j},'scramble', scram)
             %CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 1 1 2 2 2 2],'condNames',{'hit' 'cr'},'trWeights_train',[0 0 1 0 0 0], 'trWeights_test',testTrs,'num_results_iter',1)
             CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 1 1 2 2 2 2],'condNames',{'hit' 'cr'},'trWeights_train',[0 0 1 0 0 0], 'trWeights_test',testTrs,'roiName',mask{j},'scramble',scram)
         end
     end
 end

toproc = {'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te1__0__0__0__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te1__0__0__0__0__0_roirleftAngMask_tbm112514.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te1__0__0__0__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te1__0__0__0__0__0_roirleftAngMask_tbm112514.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__1__0__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__1__0__0__0_roirleftAngMask_tbm112514.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__1__0__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__1__0__0__0_roirleftAngMask_tbm112514.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__0__1_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__0__1_roirleftAngMask_tbm112514.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__0__1_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__0__1_roirleftAngMask_tbm112514.mat'

'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__1__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__1__0_roirleftAngMask_tbm112514.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__1__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__0__0__1__0_roirleftAngMask_tbm112514.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__0__1__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__0__0__1__0__0_roirleftAngMask_tbm112514.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__0__1__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__0__0__1__0__0_roirleftAngMask_tbm112514.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__1__0__0__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_exHitVexCr_trTe1_1_2_2_1_1_2_2_trWtr0__0__1__0__0__0_te0__1__0__0__0__0_roirleftAngMask_tbm112514.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__1__0__0__0__0_roiSEPT09_MVPA_MASK_resliced4mm.mat'
'conds_hitVcr_trTe1_1_1_1_2_2_2_2_trWtr0__0__1__0__0__0_te0__1__0__0__0__0_roirleftAngMask_tbm112514.mat'
}
