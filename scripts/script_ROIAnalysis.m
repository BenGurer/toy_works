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


% perform analysis and plotting seperately

% for iROI = 1:length(ROInames)
%     for iStim = 1:length(nStim)
%         saveName = ['glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim))];
%         if nStim(iStim) == 8
%             stimFreqs = stimData.stimFreqs_bin;
%             plotROIav = 0;
%             plotROIpTW = 1;
%         elseif nStim(iStim) == 32
%             stimFreqs = stimData.stimFreqs_mv;
%             plotROIav = 1;
%             plotROIpTW = 0;
%         end
%         eval([ROInames{iROI} '.roiAnalysis_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '= roiSplitGLMDataAnalysis(' ROInames{iROI} ',saveName,conditionSplitIndex,[],stimFreqs,plotROIav,plotROIpTW);'])
%     end
% end