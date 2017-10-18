function ROIdata = script_getROIdata(thisView,analysisData,analysisBaseNames,ROInames,analysisScanNum)
% ,stimInfo,plotInfo,conditionRunIndex)
% load in data
% get roi
% save data within roi to strucutre

ROIdata = struct;
q = char(39);
for iROI = 1:length(ROInames)
    eval(['ROIdata.' ROInames{iROI} ' = struct;']);
    eval(['ROIdata.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
end

% get data from ROIs
for iROI = 1:length(ROInames)
    for iAnal = 1:length(analysisBaseNames)
        %        eval([['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iScan}] '{iScan} = getROIdata_GLM(thisView,' analysisData{iScan} ',' [ 'ROIdata.' ROInames{iROI} '.roi' ] ',iScan);']);
        eval([['ROIdata.' ROInames{iROI} '.' analysisBaseNames{iAnal}] '{analysisScanNum{iAnal}} = getROIdata_GLM(thisView, analysisData{iAnal}, ROIdata.' ROInames{iROI} '.roi,analysisScanNum{iAnal});'])
        %     ANALYSIS data needs to input the analysis cell data not the strucutre scan data
        %     need to account for two analysis per scan
    end
end

%% Split GLM ROI analysis
% 
% % Perform ROI analysis
% for iROI = 1:length(ROInames)
%     for iScan = 1:length(analysisNames_scans)
%         eval([ROInames{iROI} '.roiAnalysis_' analysisBaseNames_scans{iScan} '= cal_splitAverage_roi_GLM(' ROInames{iROI} ',analysisNames_scans{iAnal},conditionRunIndex,[]);'])
%     end
% end
% % Plot ROI analysis
% for iROI = 1:length(ROInames)
%     for iPlot = 1:length(plotInfo.ROIplotList)
%         plot_splitAverage_roi_GLM(ROInames{iROI},plotInfo.ROIplotList{iPlot},stimInfo,plotInfo)
%     end
% end


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