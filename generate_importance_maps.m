for s = aubjArray
	results = res.subj{s}.penalty.nVox.weights;
	generate_importance_map(s,results,xvalIterToReport,classArgs,
end
function [impmap subj] = generate_importance_map(s,results,xvalIterToReport,class_args, numRuns)
class_args, subjId, importance_maps_dir, expt, subj, results, num_testing_sets, num?subjRuns)
results = results.iter;
expt = results.expt;
if isempty(class_args)
	class_args.train_funct_name = 'train_pLR';
end
weights = []
impmap = cell(length(expt.condNames),1);
impmap_avg = impmap;

for x = 1:length(results)
    for runNum = 1:numRuns
        switch class_args.train_funct_name
            case 'train_bp'
                weights{x}.iterations(runNum).scratchpad.net.IW{1}  = results{x}.iterations(runNum).scratchpad.net.IW{1};
            case 'train_pLR'
                weights{x}.iterations(runNum).scratchpad.net.IW{1} = results{x}.iterations(runNum).scratchpad.weights';
            case 'train_svdlr'
                weights{x}.iterations(runNum).scratchpad.net.IW{1} = results{x}.iterations(runNum).scratchpad.W';
            otherwise
                weights{x}.iterations(runNum).scratchpad.net.IW{1} = results{x}.iterations(runNum).scratchpad.w';
        end    end
end

subj = JR_extract_voxel_importance_values_kclass(subj, results,weights);
load([expt.dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
% To create vol_info.mat, run the following command in the matlab workspace:
% vol_info = spm_vol('name_of_image_file_with_same_dimensions_as_data_being_used_by_the_classifier.img'; save vol_info.mat vol_info
for condNum=1:length(expt.condNames)
    impmap{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map i
    if expt.anova_p_thresh == 1 % NO ANOVA VERSION
        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        for testSetNum = 1:num_testing_sets
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            temp{condNum}(voxel_inds)=subj.patterns{end-numRuns+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
            impmap{condNum} = impmap{condNum}+temp{condNum}; %add values cumulatively across iterations
        end
        impmap{condNum} = impmap{condNum}/numRuns*1000; %compute average and multiply by 1000 for scaling
        vol_info.fname = [importance_maps_dir '/' subjId '_' expt.condNames{condNum} '.img'];
        
    else % ANOVA VERSION
        for testSetNum = 1:num_testing_sets
            temp{condNum} = zeros(vol_info.dim); %initialize appropriately sized matrix
            voxel_inds{testSetNum} = find(subj.masks{end-numRuns+testSetNum}.mat); %get mask voxel indices
            temp{condNum}(voxel_inds{testSetNum})=subj.patterns{end-numRuns+testSetNum}.mat(:,condNum); %store impmap values at appropriate voxel indices
            impmap{condNum} = impmap{condNum}+temp{condNum}; %add values cumulatively across iterations
        end
        %sum across masks to get composite mask (where value of each voxel = number of runs for which that voxel was included)
        composite_mask = zeros(vol_info.dim);
        for maskNum = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
            composite_mask = composite_mask+subj.masks{maskNum}.mat;
        end
        voxels_to_exclude = find(composite_mask<=5);  % exclude voxels that exist for fewer than 5 of the ANOVA masks
        impmap{condNum}(voxels_to_exclude)=0;
        impmap_avg{condNum} = impmap{condNum}./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
        vol_info.fname = [importance_maps_dir '/' subjId '_' expt.condNames{condNum} '_p' num2str(expt.anova_p_thresh) '.img'];
        impmap{condNum} = impmap_avg{condNum};
    end
    spm_write_vol(vol_info,impmap{condNum});
    
    
end
end

