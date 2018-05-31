function thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,groupBase,analysisBase,overlayNumber,smoothingParams)
%% gradient reversals
% add naming to allow left and right sides for same analysis

% set what group and anaylsis to use in setup study function
%
% for iGroup = 1:length(concatenationGroup)
% for iAnalysis = 1:length(functionalAnalysis)


for iSide=1:length(subjectInfo.flatmapNames)
    % gradient reversals
    thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',subjectInfo.flatmapNames{iSide}));
    thisView = viewSet(thisView,'curgroup',groupBase);
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisBase));
    
    refreshMLRDisplay(thisView);
    params = [];
    [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlayNumber)]);
    params.combineFunction = 'gradientReversal';
    params.additionalArgs = smoothingParams;    
    params.baseSpaceInterp = 'linear';
    params.nOutputOverlays= 7;
    params.baseSpace = 1;
    params.outputName = ['gradientReversal_' Info.sides{iSide} '_' analysisBase];
    params.exportToNewGroup=1;
    [thisView,params] = combineTransformOverlays(thisView,params);
    curOverlay=viewGet(thisView,'curOverlay');
    thisView = viewSet(thisView,'overlayMin',15,curOverlay-1);
    thisView = viewSet(thisView,'overlayMax',180,curOverlay-1);
    thisView = viewSet(thisView,'overlaycolorRange',[45 180],curOverlay-1);
    thisView = viewSet(thisView,'overlayMax',75);
    thisView = viewSet(thisView,'overlaycolorRange',[0 90]);
end


%% cortical distance?
% find each voxels in an ROI cortical distance from reference
% Will probably add more