function expt = acceptUserInput(expt,  args)
%% allows user to make changes from the parameters specified in CM_mvpa_params without having to change everything every time
p = inputParser;
p.addParamValue('roiName', expt.roiName, @(x) ischar(x));
p.addParamValue('trWeights_train', expt.trWeights_train, @(x) isnumeric(x));
p.addParamValue('trWeights_test', expt.trWeights_test, @(x) isnumeric(x));
p.addParamValue('which_traintest', expt.which_traintest, @(x) isnumeric(x));
p.addParamValue('condNames', expt.condNames, @(x) iscell(x));
p.addParamValue('trnSubCondsToBal', expt.trnSubCondsToBal, @(x) iscell(x) || ischar(x));
p.addParamValue('tstSubCondsToBal', expt.tstSubCondsToBal, @(x) iscell(x) || ischar(x));
p.addParamValue('num_results_iter',expt.num_results_iter, @(x) isnumeric(x));
p.addParamValue('generate_importance_maps',expt.generate_importance_maps, @(x) (x==1 || x ==0));
p.addParamValue('xval_iters_for_imp_map',[],@(x) isnumeric(x));
p.addParamValue('scramble',expt.scramble,@(x) isnumeric(x));
p.parse(args{:});
res = p.Results;
expt = mergestructs(res, expt);
end
