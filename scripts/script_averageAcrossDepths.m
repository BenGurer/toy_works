function thisView = script_averageAcrossDepths(thisView,overlays,groupName)
    %
    %   usage: script_averageAcrossDepths(thisView,overlayNames)
    %      by: Ben Gurer
    %    date: 26/10/2017
    % purpose: average overlay across depths
    %   input: mrTools view, names of overlays, 
    %          assumes thisView is the correct group and analysis
    %

% set to flatmap group
% average all overlays

% either get all overlays or use current overlay - do i want to set overlay in here?


% use name to get overlay - use to get data
% groupName = [flatmapNames{iSide}, 'Volume'];

thisView = viewSet(thisView,'curgroup',groupName);
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
if isempty(overlays)
    analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
    overlays = 1:length(analysisData.overlays);
end


[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlays)]);
params.combineFunction='averageDepthVol';
params.combinationMode = 'Apply function to each overlay';
[thisView,params] = combineTransformOverlays(thisView,params);