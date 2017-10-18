function data = cal_splitAverage_roi_GLM(roiData,analName,conditionRunIndex,restrictIndex,conATrue)
% ,stimFreqs,plotROIav,plotROIpTW)

for iGroup = 1:length(conditionRunIndex)
    % get data from roi struct
    eval(['runA{iGroup} = roiData.' analName '{' num2str(conditionRunIndex{iGroup}(1)) '}.betas;'])
    eval(['runB{iGroup} = roiData.' analName '{' num2str(conditionRunIndex{iGroup}(2)) '}.betas;'])
    
    % Remove nans
    runA{iGroup}(isnan(runA{iGroup}(:))) = [];
    runB{iGroup}(isnan(runB{iGroup}(:))) = [];
    % restrict by index
    if ~isempty(restrictIndex)
    runA{iGroup} = runA{iGroup}(:,restrictIndex);
    runB{iGroup} = runB{iGroup}(:,restrictIndex);
    end
    
end

for iGroup = 1:length(conditionRunIndex)
    
    if size(runA{1},1) > 8
        runA{iGroup} = cal_movingAverage(runA{iGroup});
        runB{iGroup} = cal_movingAverage(runB{iGroup});
    end
    
    runA_Peak{iGroup} = cal_voxel_properties(runA{iGroup});
    runB_Peak{iGroup} = cal_voxel_properties(runB{iGroup});
    
    roi_av{iGroup} = cal_ROI_pTW_av(runA{iGroup},runB{iGroup});
    
    [condition_splitMean{iGroup}, ROI_data{iGroup}, Voxel_data{iGroup}, totalROIpCF{iGroup}] = cal_splitMean(runA{iGroup},runB{iGroup});
end

%% if condition A is true, find how pTW changes in condition B
ConB_ROIpTW = cal_ConBROIpTW_ConAVoxelIndex(Voxel_data{1}.Mean_Peak,runA{2},runB{2});

%% save data
data.roi_av = roi_av;
data.roi_pTW{1} = condition_splitMean{1};
data.roi_pCFtally{1} = totalROIpCF{1};
if conATrue == 1
data.roi_pTW{2} = ConB_ROIpTW;
data.roi_pCFtally{2} = {nan};
else
data.roi_pTW{2} = condition_splitMean{2};
data.roi_pCFtally{2} = totalROIpCF{2};
end
