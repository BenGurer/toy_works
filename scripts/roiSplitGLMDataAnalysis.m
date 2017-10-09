function data = roiSplitGLMDataAnalysis(roiData,analName,conIndex,restrictIndex)

for iGroup = 1:length(conIndex)
    % get data from roi struct
    eval(['runA{iGroup} = roiData.' analName '{' num2str(conIndex{iGroup}(1)) '}'])
    eval(['runB{iGroup} = roiData.' analName '{' num2str(conIndex{iGroup}(2)) '}'])
    
    % Remove nans
    runA{iGroup}(isnan(runA{iGroup})) = [];
    runB{iGroup}(isnan(runB{iGroup})) = [];
    % restrict by index
    runA{iGroup} = runA{iGroup}(:,restrictIndex);
    runB{iGroup} = runB{iGroup}(:,restrictIndex);
    
end

for iGroup = 1:length(conIndex)
    
    if length(conditionA) > 8
    runA{iGroup} = cal_movingAverage(runA{iGroup});
    runB{iGroup} = cal_movingAverage(runB{iGroup});
    end
    
    runA_Peak_mv{iGroup} = cal_voxel_properties(runA{iGroup});
    runB_Peak_mv{iGroup} = cal_voxel_properties(runB{iGroup});
    
    roi_av_mv{iGroup} = cal_ROI_pTW_av(runA{iGroup},runB{iGroup});
    
    [condition_splitMean_bin{iGroup}, ROI_data_bin{iGroup}, Voxel_data_bin{iGroup}, totalROIpCF_bin{iGroup}] = cal_splitMean(runA{iGroup},runB{iGroup});
end

%% compare ROI average Beta weights
plot_compareConditions_ROIav(roi_av_bin{1},roi_av_bin{2},stimFreqs_bin)

%% Compare average voxel tuning
% Compare ROI pTW
plot_compareConditions_pTW(condition_splitMean_mv{1},condition_splitMean_mv{2},4);

% Compare ROI pTW assuming condition A is TRUE
ConBROIpTW_bin = cal_ConBROIpTW_ConAVoxelIndex(Voxel_data_bin{1}.Mean_Peak,runA_8bins{2},runB_8bins{2});
plot_compareConditions_pTW(condition_splitMean_bin{1}, ConBROIpTW_bin,2,stimFreqs_bin);

end