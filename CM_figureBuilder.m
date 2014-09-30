function CM_figureBuilder()
trs_to_use = 1:6;
res_dir = '/Volumes/Tyler_Drive1/CM_mvpa_results';
%first, let's look at ex>ex at each tr
res_arr = {'conds_hitVcr_trTe1_2_3_4_0_0_0_0_trW' 'conds_hitVcr_trTe0_0_0_0_1_2_3_4_trW'};
tr_str_arr = {'1__0__0__0__0__0' '0__1__0__0__0__0' '0__0__1__0__0__0'   '0__0__0__1__0__0' '0__0__0__0__1__0' '0__0__0__0__0__1'}; 
mask_str = 'roiSEPT09_MVPA*';

colors = {'g','r'};

f =figure;

for r = 1:length(res_arr)
    res_str = res_arr{r};
    for tr = trs_to_use
        fname = dir([res_dir '/' res_str tr_str_arr{tr} '_' mask_str 'mat']);
        res= load(fullfile(res_dir, fname.name));
        
        res = res.res;
        
        proc_res{tr} = CM_mvpaPostProc(res, [], 'ret', 0, '',[]);
        auc(tr,:) = proc_res{tr}.auc
        
    end
    
    mean_auc = nanmean(auc,2);
    ste_auc = nansem(auc,2);
    
    set(0, 'currentfigure', f);  
    hold on;
    errorbar3(trs_to_use,mean_auc',ste_auc',1,colors{r});
    hold on;
    plot(trs_to_use,mean_auc',colors{r}, 'LineWidth',1.5);
end

axis square;
end



