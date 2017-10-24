function data = script_covertData2FlatmapSpace(groupName,analysisName,overlays,flatmapNames)

for iSide=1:2
  %first copy overlays to flat volume
  thisView = viewSet(thisView,'curgroup',groupName);
  thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapName{iSide}));
  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlays)]);
  params.combineFunction='User Defined';
  params.customCombineFunction = 'plus'; %add 0
  params.combinationMode = 'Apply function to each overlay';
  params.additionalArgs = '0';
  params.outputName=' ';
  params.baseSpace = 1; %export result to base space (flat map)
  params.exportToNewGroup=1; %export to volume in a new group
  [thisView,params] = combineTransformOverlays(thisView,params);
%   curOverlay=viewGet(thisView,'curOverlay');
%   thisView = viewSet(thisView,'overlayColorRange',[0 5],curOverlay-[0 1 3]);
%   thisView = viewSet(thisView,'overlayColorRange',[0 8],curOverlay-2);
%   thisView = viewSet(thisView,'curSlice',6);
%   curOverlay=viewGet(thisView,'curOverlay');
%   thisView = viewSet(thisView,'overlaycolorrange',[0 8],curOverlay-2);
%   thisView = viewSet(thisView,'overlaycolorrange',[0 5],curOverlay-1);
%   thisView = viewSet(thisView,'overlaycolorrange',[0 5],curOverlay);
%   thisView = viewSet(thisView,'overlayMin',1,curOverlay-3);
%   thisView = viewSet(thisView,'alphaOverlay',curOverlay-3,curOverlay-(0:2));
%   thisView = viewSet(thisView,'alphaOverlayExponent',0,curOverlay-(0:2));
end