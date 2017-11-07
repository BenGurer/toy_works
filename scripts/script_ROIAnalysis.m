function roiData = script_ROIAnalysis(roiData,analysisNames,Info,stimInfo,plotInfo,conditionRunIndex,analysisScanNum,dataType,ROInamesSelect)

%% Split GLM ROI analysis
% Perform ROI analysis
if isempty(ROInamesSelect)
    ROInames = Info.ROInames;
else
    if iscell(ROInamesSelect)
        ROInames = ROInamesSelect;
    else
        ROInames = {ROInamesSelect};
    end
end
%% perform ROI analysis on GLM data
for iROI = 1:length(ROInames)
    for iAnal = 1:length(analysisNames)
        eval(['roiData.' ROInames{iROI} '.roiAnalysis_' analysisNames{iAnal} '= cal_splitAverage_roi_GLM(roiData.' ROInames{iROI} ',analysisNames{iAnal},conditionRunIndex,[],Info.ConATrue,dataType);'])
    end
end
%% calculate ratio between dB senstation level and BOLD fMRI activity




%% Plot ROI analysis
for iROI = 1:length(ROInames)
    for iPlot = 1:length(plotInfo.ROIplotList)
        plot_splitAverage_roi_GLM(eval(['roiData.' ROInames{iROI}]),plotInfo.ROIplotList{iPlot},stimInfo,plotInfo,ROInames{iROI});
    end
end

