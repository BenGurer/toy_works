function data = script_averageAcrossDepths



% set to flatmap group
% average all overlays

% use name to get overlay - use to get data
a = viewGet(thisView,'Overlay','Scan 1 -  (Tone 100Hz,0)');
a = viewGet(thisView,'Overlay',['Scan 1 -  (Tone ', * 'Hz,0)']);

groupName = [flatmapNames{iSide}, 'Volume'];

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