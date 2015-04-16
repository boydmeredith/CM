function [] = auc_wrap()
scram_run1 = load('scrambled_conds_exHitVexCr_trTe1_2_2_2_1_1_1_1_trWtr0___________0________0.33________0.34________0.33___________0_te0___________0________0.33________0.34________0.33___________0_roiSEPT09_MVPA_MASK_resliced4mm.mat');
unscram_run1 = load('conds_exHitVexCr_trTe1_2_2_2_1_1_1_1_trWtr0___________0________0.33________0.34________0.33___________0_te0___________0________0.33________0.34________0.33___________0_roiSEPT09_MVPA_MASK_resliced4mm.mat');
runs_1 = 2:4;
unscram_run4 = load('conds_exHitVexCr_trTe2_2_2_1_1_1_1_1_trWtr0___________0________0.33________0.34________0.33___________0_te0___________0________0.33________0.34________0.33___________0_roiSEPT09_MVPA_MASK_resliced4mm.mat');
scram_run4 = load('scrambled_conds_exHitVexCr_trTe2_2_2_1_1_1_1_1_trWtr0___________0________0.33________0.34________0.33___________0_te0___________0________0.33________0.34________0.33___________0_roiSEPT09_MVPA_MASK_resliced4mm.mat');
runs_4 = 1:3;

[auc_1_un auc_1_scram] = auc_by_run(unscram_run1, scram_run1, runs_1)
[auc_4_un auc_4_scram] = auc_by_run(unscram_run4, scram_run4, runs_4)
keyboard

function [unscram_auc scram_auc] = auc_by_run(unscram, scram, runs)

subNos = [1 3:10 12:26];
iters = 1:10;
unscram_auc = nan(size(iters),size(subNos),size(runs));
scram_auc = nan(size(iters),size(subNos),size(runs));
for s = subNos
for i = iters
c_s = scram.res.subj{s}.penalty.nVox.weights.condensed_runs
c = unscram.res.subj{s}.penalty.nVox.weights.condensed_runs
for r = runs
t_idx_un = unscram.res.subj{s}.penalty.nVox.weights.iter{i}.iterations(2).test_idx
t_idx_un = find(ismember(t_idx_un, find(c_s==r)))
t_idx_s = scram.res.subj{s}.penalty.nVox.weights.iter{i}.iterations(2).test_idx
t_idx_s = find(ismember(t_idx_s, find(c==r)))

des = unscram.res.subj{s}.penalty.nVox.weights.iter{i}.iterations(2).perfmet.desireds(t_idx_un)
des_s = scram.res.subj{s}.penalty.nVox.weights.iter{i}.iterations(2).perfmet.desireds(t_idx_s)
acts_s = scram.res.subj{s}.penalty.nVox.weights.iter{i}.iterations(2).acts(:,t_idx_s)
acts = unscram.res.subj{s}.penalty.nVox.weights.iter{i}.iterations(2).acts(:, t_idx_un)
scram_auc(i,s,r) = calculate_auc(acts_s, des_s, 0)
unscram_auc(i,s,r) = calculate_auc(acts, des, 0)
end
end
end
