function analysisData = getGroupAnalysisData(thisView,analysisName,iScan)
% if iScan ~= 0
% else
%     iScan = 1;
% end
% get analysis data for specified group
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
% get analysis data for specified analysis
analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
% glmData = analysisData.d{iScan};
% params= analysisData.params;
end