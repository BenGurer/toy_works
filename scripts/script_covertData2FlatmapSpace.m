function data = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapNames)

if ~isempty(iScan)
    thisView = viewSet(thisView,'curgroup',groupName,['curScan=' mat2str(iScan)]);
else
    thisView = viewSet(thisView,'curgroup',groupName);
end
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
if isempty(overlays)
    analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
    overlays = 1:length(analysisData.overlays);
end

data.conditionNames = analysisData.params.EVnames;
for i = 1:length(analysisData.params.EVnames)
    if ~isempty(iScan)
        data.overlayConditionNames{i} = ['Scan ', mat2str(iScan), ' -  (' analysisData.params.EVnames{i} ',0)'];
    else
        data.overlayConditionNames{i} = [groupName ' (' analysisData.params.EVnames{i} ',0)'];
    end
end

a = viewGet(thisView,'Overlay',data.overlayConditionNames{1})

for iSide=1:2
    % copy overlays to flat volume
    thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',flatmapNames{iSide}));
    [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlays)]);
    params.combineFunction='User Defined';
    params.customCombineFunction = 'plus'; %add 0
    params.combinationMode = 'Apply function to each overlay';
    params.additionalArgs = '0';
    if ~isempty(iScan)
        params.outputName=' ';
    else
        params.outputName= [groupName ' '];
    end
    params.baseSpace = 1; %export result to base space (flat map)
    params.exportToNewGroup=1; %export to volume in a new group
    [thisView,params] = combineTransformOverlays(thisView,params);
end