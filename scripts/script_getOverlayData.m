function overlayData = script_getOverlayData(thisView,groupName,analysisName,conditionNames,iScan)
    %
    %   usage: script_getOverlayData(thisView,groupName)
    %      by: Ben Gurer
    %    date: 26/10/2017
    % purpose: Get overlay data from  group
    %   input: mrTools view, names of overlays

    if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
        thisView = viewSet(thisView,'curgroup',groupName);
    end
    
    if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
    end
    
% thisView = viewSet(thisView,'curgroup',groupName);
% thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
% if isempty(overlays)
%     analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
%     overlays = 1:length(analysisData.overlays);
% end
% if ~isempty(iScan)
%     for iCon = 1:length(conditionNames)
%         overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' analysisName ' (' conditionNames{iCon} ',0))'];
%     end
% overlayData = get_overlayData(thisView,overlayNames);
% else
    overlayData = get_overlayData(thisView,conditionNames);
% end
end