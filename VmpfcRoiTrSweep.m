masksToClassify= {'rleftIfgAnat.nii', 'rleftBa44_45_47anat.nii'};
%classify late trs for all rois
for scram = 0:1
    for j = 1:length(masksToClassify)
     for i = [1:6]
             testTrs = zeros(1,6);
             testTrs(i) = 1
        CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 2 2 1 1 2 2 ],'condNames',{'exHit' 'exCr'},'scramble', scram,'trWeights_train',[0 0 1 0 0 0], 'trWeights_test',testTrs,'roiName',masksToClassify{j})
        CM_run_mvpa_v4_tbm([],'','which_traintest',[1 1 1 1 2 2 2 2],'condNames',{'hit' 'cr'},'scramble',scram,'trWeights_train',[0 0 1 0 0 0], 'trWeights_test',testTrs,'roiName',masksToClassify{j})
	end
	end
end

%process it all
cd('~/cm/Data/Functional/mvpa_results/')
%all late rois
CM_batchMvpaPostProc({'*ex*ex*leftIfgAnat*mat' '*hitVcr*leftIfgAnat*mat' '*ex*ex*Ba44_45*mat' '*hitVcr*Ba44_45*mat'}, 'pfc_rois_trsweep',pwd,0, {[],2, [], 2})
%trbytr whole brain and motor
CM_batchMvpaPostProc({'*ex*ex*1_2_3_4_1_2_3_4*locer*mat' '*ex*ex*1_2_3_4_1_2_3_4*SEPT*mat' '*hitVcr*1_1_1_1_2_2_2_2*locer*mat' '*hitVcr*1_1_1_1_2_2_2_2*SEPT*mat'}, 'wb_motor_trbytr',pwd,0, {[],[],2, 2})
