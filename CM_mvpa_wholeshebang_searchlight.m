AA_subjs = {'CM004', 'CM008', 'CM010', 'CM013', 'CM019', 'CM020', 'CM021', 'CM025'};
all_subjs = {'CM001','CM003','CM004','CM005','CM006','CM007','CM008','CM009','CM010','CM012','CM013','CM014','CM015','CM016','CM017','CM018','CM019','CM020','CM021','CM022','CM023','CM024','CM025','CM026'};
highEXtoEXClassPerfSubs = {'CM001' ,   'CM005' ,   'CM006'  ,  'CM009' ,   'CM010'   , 'CM012'   , 'CM014' , 'CM015' ,   'CM019' ,   'CM022'    ,'CM023' ,   'CM025'};

highClassPerfSubIdx = [1     4     5     8     9    10    12    13    17    20    21    23];

loTaskDiffSubs = [2     3    13    14    15    16    17    18    20    21    23    24];

%  CM_run_mvpa_v3_tbm_kclass(all_subjs, {'EX_hits','EX_CRs', 'CM_hits', 'CM_CRs'}, 0, 10, 7, 345, 'rFaceGrScene05k5ANDBLFusiformAnat');
%  CM_run_mvpa_v3_tbm_kclass(all_subjs, {'EX_hits','EX_misses', 'EX_CRs', 'EX_FAs'}, 0, 10, 1, 345, 'rleft_VLPFC_anat');
% CM_run_mvpa_v3_tbm_kclass(all_subjs, {'EX_hits','EX_CRs', 'CM_hits', 'CM_CRs'}, 0, 10, 7, 345, 'rAnG_anat');
% CM_run_mvpa_v3_tbm_kclass(all_subjs, {'EX_hits','EX_misses', 'EX_CRs', 'EX_FAs'}, 0, 10, 1, 345, 'rAnG_anat');

% CM_run_mvpa_v3_tbm_kclass(all_subjs, {'EX_hits','EX_FAs', 'CM_hits', 'CM_FAs'}, 0, 10, 7, 345, '');
% CM_run_mvpa_v3_tbm_kclass(all_subjs, {'EX_hits','EX_FAs', 'CM_hits', 'CM_FAs'}, 0, 10, 7, 345, '');

% CM_run_mvpa_v3_tbm(all_subjs,'EX_hits','EX_FAs',0,10,3,3)
%CM_run_mvpa_searchlight_batch(all_subjs, 'EX_hits','EX_CRs',1,10,10,'ver1');
% rad = 1;
% for i =2:5
%     nickname = sprintf('ver%i',i)
%     CM_run_mvpa_searchlight_batch(all_subjs, 'EX_hits','EX_CRs',rad,10,1,nickname);
%     CM_run_mvpa_searchlight_batch(all_subjs, 'CM_hits','CM_CRs',rad,10,2,nickname);
%     CM_run_mvpa_searchlight_batch(all_subjs, 'EX_and_CM_hits','EX_and_CM_CRs',rad,10,10,nickname);
% end
% for rad = [ 7];
%     for i = 1
%         nickname = sprintf('ver%i',i)
%         CM_run_mvpa_searchlight_batch(thistimesubs, 'EX_hits','EX_CRs',rad,10,1,nickname);
%         %CM_run_mvpa_searchlight_batch(all_subjs, 'CM_hits','CM_CRs',rad,10,2,nickname);
%         %     CM_run_mvpa_searchlight_batch(all_subjs, 'EX_and_CM_hits','EX_and_CM_CRs',rad,10,10,nickname);
%     end
% end
for rad = [8 ]
    CM_run_mvpa_searchlight_batch(all_subjs, {'EX_hits','EX_CRs'},rad,10,1,'');
    CM_run_mvpa_searchlight_batch(all_subjs, {'CM_hits','CM_CRs'},rad,10,2,'');
end
    