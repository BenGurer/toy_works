function data = cal_average_roi_GLM(roiData, groupNames, analName,restrictIndex,conATrue,dataType,stimInfo)
% ,stimFreqs,plotROIav,plotROIpTW)
q = char(39);
for iGroup = 1:length(groupNames)
    % get data from roi struct
    % convert cell to mat and save in betas before this
    % use input arguement to set or set before to be run a and b
    %     switch dataType
    %         case 'GLM'
    %     eval(['runA{iGroup} = roiData.' analName '{' num2str(groupNames{iGroup}) '}.betas;'])
    %     eval(['runB{iGroup} = roiData.' analName '{' num2str(groupNames{iGroup}) '}.betas;'])
    %         case 'overlays'
%     eval(['raw{iGroup} = cell2mat(roiData.' groupNames{iGroup} '.' roiName '.' analName q ');'])
    %     eval(['conB = cell2mat(roiData.' groupNames{iGroup} '.' roiName '.' analName ');'])
    %     end
   
    eval(['raw{iGroup} = cell2mat(roiData.' groupNames{iGroup} '.' analName q ');']) 
    
    % Remove nans
    raw{iGroup}(isnan(raw{iGroup}(:))) = [];
    % restrict by index
    if ~isempty(restrictIndex)
        raw{iGroup} = raw{iGroup}(:,restrictIndex);
    end
    
    raw_mv{iGroup} = cal_movingAverage(raw{iGroup});
    
    roi_av{iGroup} = mean(raw_mv{iGroup},2);
    
    raw_Peak{iGroup} = cal_voxel_properties(raw{iGroup});
    
end

% %% if condition A is true, find how pTW changes in condition B
% ConB_ROIpTW = cal_ConBROIpTW_ConAVoxelIndex(raw_Peak{1},runA{2},runB{2});

%% save data
data.roi_av = roi_av;
% data.roi_pTW{1} = condition_splitMean{1};
data.roi_pCFtally = raw_Peak;
data.roi_av_ratio = roi_av{2}./roi_av{1};
% if conATrue == 1
%     data.roi_pTW{2} = ConB_ROIpTW;
%     data.roi_pCFtally{2} = {nan};
% else
%     data.roi_pTW{2} = condition_splitMean{2};
%     data.roi_pCFtally{2} = totalROIpCF{2};
% end

if length(roi_av{1}) == stimInfo.sizes(1)
    stimulusLevels = stimInfo.stimLevel_SL_bin;
elseif length(roi_av{1}) == stimInfo.sizes(2)
    stimulusLevels = stimInfo.stimLevel_SL_mv;
else 
    stimulusLevels = stimInfo.stimLevel_SL;
end
baseLevel_dB = 50;
[ data.ratio2Plot, data.level2Plot, data.fit , data.error2plot ] =  cal_dbSLvsBetaWeight(data.roi_av_ratio,stimulusLevels,baseLevel_dB);
