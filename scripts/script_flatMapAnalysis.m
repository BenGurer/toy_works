function thisView = script_flatMapAnalysis(thisView,Info,subjectInfo)
%% gradient reversals
% add naming to allow left and right sides for same analysis

% set what group and anaylsis to use in setup study function
%
% for iGroup = 1:length(concatenationGroup)
% for iAnalysis = 1:length(functionalAnalysis)
for iSide=1:length(Info.sides)
    % gradient reversals
    thisView = viewSet(thisView,'curgroup',Info.gradReversalInfo.groupBase);
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',Info.gradReversalInfo.analysisBase));
    %     thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['x' subjectInfo.freeSurferName '_' Info.sides{iSide} '_WM_Flat_' subjectInfo.flatmapNames{iSide}]));
    thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',subjectInfo.flatmapNames{iSide}));
    
    refreshMLRDisplay(thisView);
    [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(Info.gradReversalInfo.overlayBase)]);
    params.combineFunction = 'gradientReversal';
    params.additionalArgs = '[18 18 21]';
    params.baseSpaceInterp = 'linear';
    params.nOutputOverlays= 7;
    params.baseSpace = 1;
    params.outputName = ['gradientReversal_' Info.sides{iSide}];
    params.exportToNewGroup=1;
    [thisView,params] = combineTransformOverlays(thisView,params);
    curOverlay=viewGet(thisView,'curOverlay');
    thisView = viewSet(thisView,'overlayMin',15,curOverlay-1);
    thisView = viewSet(thisView,'overlayMax',180,curOverlay-1);
    thisView = viewSet(thisView,'overlaycolorRange',[45 180],curOverlay-1);
    thisView = viewSet(thisView,'overlayMax',75);
    thisView = viewSet(thisView,'overlaycolorRange',[0 90]);
end
% end
% end

%% cortical distance?
% find each voxels in an ROI cortical distance from reference
% Will probably add more