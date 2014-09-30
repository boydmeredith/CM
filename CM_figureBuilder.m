function CM_figureBuilder()
trs_to_use = 1:6;
res_dir = '/Volumes/Tyler_Drive1/CM_mvpa_results';
%first, let's look at ex>ex at each tr
res_str = 'conds_hitVcr_trTe1_2_3_4_0_0_0_0_trW';
tr_str = {'0__0__0__0__0__1' '0__0__0__0__1__0' '0__0__0__1__0__0' '0__0__1__0__0__0' '0__1__0__0__0__0' '1__0__0__0__0__0'}
mask_str = 'roiSEPT09_MVPA*';

colors = {};

figure;


for tr = trs_to_use
    fname = dir([res_dir '/' res_str tr_str{tr} '_' mask_str 'mat']);
    res= load(fullfile(res_dir, fname.name)); 
    
    res = res.res;
    
    proc_res{tr} = CM_mvpaPostProc(res, [], 'ret', '',[]);
    aucs(tr,:) = proc_res{tr}.auc
        
end

plot(trs_to_use,nanmean(auc,2));

errorbar3(trs_to_use,nanmean(auc,2),ste(auc,2),1,[.2 .6 1]);

end


    
    