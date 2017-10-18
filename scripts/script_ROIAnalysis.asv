function roiData = script_ROIAnalysis(roiData,analysisNames,Info,stimInfo,plotInfo,conditionRunIndex,analysisScanNum)

%% Split GLM ROI analysis
% Perform ROI analysis
for iROI = 1:length(Info.ROInames)
    for iAnal = 1:length(analysisNames)
        eval(['roiData.' Info.ROInames{iROI} '.roiAnalysis_' analysisNames{iAnal} '= cal_splitAverage_roi_GLM(roiData.' Info.ROInames{iROI} ',analysisNames{iAnal},conditionRunIndex,[],Info.ConATrue);'])
    end
end
% Plot ROI analysis
for iROI = 1:length(Info.ROInames)
    for iPlot = 1:length(plotInfo.ROIplotList)
        plot_splitAverage_roi_GLM(eval(['roiData.' Info.ROInames{iROI}]),plotInfo.ROIplotList{iPlot},stimInfo,plotInfo);
    end
end
