function data = cal_splitAverage_roi_GLM(roiData,analName,conditionRunIndex,restrictIndex,conATrue,dataType,stimInfo)
% ,stimFreqs,plotROIav,plotROIpTW)
q = char(39);
for iGroup = 1:length(conditionRunIndex)
    % get data from roi struct
    % convert cell to mat and save in betas before this
    % use input arguement to set or set before to be run a and b
    switch dataType
        case 'GLM'
    eval(['runA{iGroup} = roiData.' analName '{' num2str(conditionRunIndex{iGroup}(1)) '}.betas;'])
    eval(['runB{iGroup} = roiData.' analName '{' num2str(conditionRunIndex{iGroup}(2)) '}.betas;'])
        case 'overlays'            
    eval(['runA{iGroup} = cell2mat(roiData.' analName '{' num2str(conditionRunIndex{iGroup}(1)) '}' q ');'])
    eval(['runB{iGroup} = cell2mat(roiData.' analName '{' num2str(conditionRunIndex{iGroup}(2)) '}' q ');'])            
    end
    
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
data.roi_av_ratio = roi_av{2}./roi_av{1};
if conATrue == 1
data.roi_pTW{2} = ConB_ROIpTW;
data.roi_pCFtally{2} = {nan};
else
data.roi_pTW{2} = condition_splitMean{2};
data.roi_pCFtally{2} = totalROIpCF{2};
end


if length(roi_av{1}) == stimInfo.sizes(1)
    stimulusLevels = stimInfo.stimLevel_SL_bin;
elseif length(roi_av{1}) == stimInfo.sizes(2)
    stimulusLevels = stimInfo.stimLevel_SL_mv;
else 
    stimulusLevels = stimInfo.stimLevel_SL;
end
baseLevel_dB = 50;
[ data.ratio2Plot, data.level2Plot, data.fit , data.error2plot ] =  cal_dbSLvsBetaWeight(data.roi_av_ratio,stimulusLevels,baseLevel_dB);
