function data = roiSplitGLMDataAnalysis(roiData,analName,conIndex,restrictIndex)
% ,stimFreqs,plotROIav,plotROIpTW)

for iGroup = 1:length(conIndex)
    % get data from roi struct
    eval(['runA{iGroup} = roiData.' analName '{' num2str(conIndex{iGroup}(1)) '}.betas;'])
    eval(['runB{iGroup} = roiData.' analName '{' num2str(conIndex{iGroup}(2)) '}.betas;'])
    
    % Remove nans
    runA{iGroup}(isnan(runA{iGroup}(:))) = [];
    runB{iGroup}(isnan(runB{iGroup}(:))) = [];
    % restrict by index
    if ~isempty(restrictIndex)
    runA{iGroup} = runA{iGroup}(:,restrictIndex);
    runB{iGroup} = runB{iGroup}(:,restrictIndex);
    end
    
end

if size(runA{1},1) <= 8
    nCols = 2;
else
    nCols = 4;
end

for iGroup = 1:length(conIndex)
    
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

% if plotROIav == 1
% %% compare ROI average Beta weights
% plot_compareConditions_ROIav(roi_av{1},roi_av{2},stimFreqs)
% 
% end
% 
% if plotROIpTW == 1
% %% Compare average voxel tuning
% % Compare ROI pTW
% % plot_compareConditions_pTW(condition_splitMean{1},condition_splitMean{2},nCols,stimFreqs);
% 
% % Compare ROI pTW assuming condition A is TRUE
% plot_compareConditions_pTW(condition_splitMean{1}, ConB_ROIpTW,nCols,stimFreqs);
% 
% end

%% save data
data.roi_av = roi_av;
data.roi_pTW{1} = condition_splitMean{1};
data.roi_pTW{2} = ConB_ROIpTW;
end