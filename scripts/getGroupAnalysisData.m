function analysisData = getGroupAnalysisData(thisView,analysisName)
    %
    %   usage: getGroupAnalysisData(thisView,analysisName,iScan)
    %      by: Ben Gurer
    %    date: 05/10/2017
    % purpose: get analysis data from mrTools view
    %   input: mrTools view, analysis name, ROI, scan number (if defined)
    %
% get analysis data for specified group
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
% get analysis data for specified analysis
analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
end