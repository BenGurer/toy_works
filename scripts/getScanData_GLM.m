function data = getScanData_GLM(thisView,analysisNames_Scans,analysisNames_Groups,groupNames)
%% Get data from individual scans
thisView = viewSet(thisView,'curGroup','MotionComp');
% for iStim = 1:length(nStim)
%     for iScan = 1:nScans
for iScan = 1:length(analysisNames_Scans)
    %                eval(['GLManalysisData_nCons_' mat2str(nStim(iStim)) '{iScan} = getGroupAnalysisData(thisView,analysisName);'])
    
%     thisView = viewSet(thisView,'curGroup','MotionComp',['curScan=' mat2str(iScan)]);
    eval(['data.scan_GLMdata{iScan} = getGroupAnalysisData(thisView,analysisNames_Scans{iScan});'])

end
% end

% get data from concatenated groups
for iGroup = 1:length(groupNames)
    thisView = viewSet(thisView,'curGroup',groupNames{iGroup});
    for iAnal = 1:length(analysisNames_Groups)
        eval(['data.' groupNames{iGroup} '_GLMdata{iAnal} = getGroupAnalysisData(thisView,analysisNames_Groups{iAnal});'])
    end
    %
    %     ROIEstimatesData_Concat{iGroup} = getGroupROIEstimates(thisView,pacROI,concatenationGroupNames{iGroup},['glm_' hrfModel{2}],0);
end
%%
% get experimental data