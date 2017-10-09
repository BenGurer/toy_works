function data = roiSplitGLMDataAnalysis(roiData)

runA_32cons = [];
runB_32cons = [];
runA_8bins = [];
runB_8bins = [];
runA_32cons_restricted = [];
runB_32cons_restricted = [];
runA_8bins_restricted = [];
runB_8bins_restricted = [];

runA_mv = [];
runB_mv = [];
for iGroup = 1:length(concatenationGroupNames)
    runA_32cons{iGroup} = ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(1)}.betas;
    runB_32cons{iGroup} = ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(2)}.betas;
    
    runA_32cons{iGroup}(isnan(runA_32cons{iGroup})) = [];
    runB_32cons{iGroup}(isnan(runB_32cons{iGroup})) = [];
    
    runA_8bins{iGroup} = ROIEstimatesData_SplitRuns_8bins{conditionScans{iGroup}(1)}.betas;
    runB_8bins{iGroup} = ROIEstimatesData_SplitRuns_8bins{conditionScans{iGroup}(2)}.betas;
    
    runA_32cons_restricted{iGroup} = runA_32cons{iGroup}(:,restrictIndex);
    runB_32cons_restricted{iGroup} = runB_32cons{iGroup}(:,restrictIndex);
    runA_8bins_restricted{iGroup} = runA_8bins{iGroup}(:,restrictIndex);
    runB_8bins_restricted{iGroup} = runB_8bins{iGroup}(:,restrictIndex);
    
    runA_mv{iGroup} = cal_movingAverage(runA_32cons_restricted{iGroup});
    runB_mv{iGroup} = cal_movingAverage(runB_32cons_restricted{iGroup});
    
    runA_Peak_mv{iGroup} = cal_voxel_properties(runA_mv{iGroup});
    runB_Peak_mv{iGroup} = cal_voxel_properties(runB_mv{iGroup});
    
    roi_av{iGroup} = cal_ROI_pTW_av(runA_32cons_restricted{iGroup},runB_32cons_restricted{iGroup});
    roi_av_mv{iGroup} = cal_ROI_pTW_av(runA_mv{iGroup},runB_mv{iGroup});
    roi_av_bin{iGroup} = cal_ROI_pTW_av(runA_8bins_restricted{iGroup},runB_8bins_restricted{iGroup});
    
    % [condition_splitMean{iGroup}, ROI_data{iGroup}, Voxel_data{iGroup}, totalROIpCF{iGroup}] = cal_splitMean(runA_32cons_restricted{iGroup},runB_32cons_restricted{iGroup});
    [condition_splitMean_mv{iGroup}, ROI_data_mv{iGroup}, Voxel_data_mv{iGroup}, totalROIpCF_mv{iGroup}] = cal_splitMean(runA_mv{iGroup},runB_mv{iGroup});
    [condition_splitMean_bin{iGroup}, ROI_data_bin{iGroup}, Voxel_data_bin{iGroup}, totalROIpCF_bin{iGroup}] = cal_splitMean(runA_8bins_restricted{iGroup},runB_8bins_restricted{iGroup});
end

%% compare ROI average Beta weights

%% compare 



end