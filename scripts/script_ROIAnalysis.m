function data = script_ROIAnalysis(thisView)

%% Save data in structures labelled by roi
%% ROI names
% Create ROIs with these names:
% RightAC
% RightPosAC
% RightAntAc
% LeftAC
% LeftPosAC
% LeftAntAC
% AC
ROInames = {'RightAC','RightPosAC','RightAntAC','LeftAC','LeftPosAC','LeftAntAC','AC'};

ROInames = {'RightAC','RightPosAC'};

q = char(39);
for iROI = 1:length(ROInames)
eval([ROInames{iROI} ' = struct;']);
eval([ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
end
% function to create figure based on paper size - poster or paper
% send figure handle to plotting functions - can be used with subplot and plot
% create strucutre for rois or containing rois

% following analysis - save roi results underone structure of subject ID
% plots I want
% ROI average beta weight
% ROI TW conA v B

%% SHOULD SET ANALYSIS NAME ONCE AND PASS IT ON
% loop - create analysis names
% loop length of analysis names
% needs to name based on group - one structure

%% Get data from individual scans
thisView = viewSet(thisView,'curGroup','MotionComp');
for iStim = 1:length(nStim)
    for iScan = 1:nScans
        analysisName = ['glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '_Scan_' mat2str(iScan)];
        eval(['GLManalysisData_nCons_' mat2str(nStim(iStim)) '{iScan} = getGroupAnalysisData(thisView,analysisName);'])
    end
end
% get data from ROIs
for iROI = 1:length(ROInames)
    for iScan = 1:nScans
        for iStim = 1:length(nStim)
            eval([ROInames{iROI} '.glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '{iScan} = getROIdata(thisView,GLManalysisData_nCons_' mat2str(nStim(iStim)) '{iScan},' ROInames{iROI} '.roi,iScan);']);
        end
    end
end

% now get concat group data
thisView = viewSet(thisView,'curGroup',groupName);

% get datga from concatenated groups
for iGroup = 1:length(concatenationGroupNames)
    
thisView = viewSet(thisView,'curGroup',concatenationGroupNames{iGroup});

    ROIEstimatesData_Concat{iGroup} = getGroupROIEstimates(thisView,pacROI,concatenationGroupNames{iGroup},['glm_' hrfModel{2}],0);
end
%%
% get experimental data

%% now pass data to function to plot
data = roiDataAnalysis(ROIdata);
% loop over ROIs
