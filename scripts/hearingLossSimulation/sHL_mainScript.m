%% sHL_mainScript - scripted analysis for Simulated Hearing Loss (sHL) study 2017
% simulate hearing loss

%% run this function when returning for analysis?
iSub = 1;
thisView = getMLRView;
[stimInfo, glmInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

% define and get subject info
% either cell for each subject or return single subject info - prob the later
subjectInfo = get_SubjectInfo_sHL(iSub);
% Subject ID
% flatmap names

% setup study direction (once)
% create folders needed
sHL_createStudyDirectory

% import, convert and move data
sHL_organiseData(studyDirectory,subject)
% per subject

preprocessData = sHL_preprocess;
% distortion correct
% linear alignment
% non-linear alignment

[thisView, concatedate] = script_setupmrLoadRet(thisView,groupNames);
% initiate mrLoadRet
% motion corerection
% group data

[thisView, glmData] = script_glmAnalysis(thisView);
% HRF = double gamma and box car
% All stims and 8 bins
% don't need to perform weighted mean on individual runs

thisView = script_importAnatomy(thisView);
% load in:
% reference EPI
% High resolution in-plane T2*
% surfaces
% create flatmaps
% MAKE ORIGINAL FLAT MAPS (LEFT AND RIGHT) USING MAKEFLAT AND NAME THEM
% [freeSurferName{iSubj} '_left_Flat.off'] AND [freeSurferName{iSubj} '_right_Flat.off']

[thisView, pRFdata] = script_pRFAnalysis(thisView);

thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);

% get condition names
conditionNames = cell(1,length(glmInfo.nStim));
for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
    analysisName = glmInfo.analysisNames_Scans{iAnal};
    conditionNames{iAnal} = get_analysisConditionNames(thisView,analysisName,'MotionComp',1);
end

%% Create subject data storage structure

%% Get analysis data from scans
% pass on to ROI analysis to restrict

%% get data from scans and groups
% save so don't need to load again
% change name to get_analysisData_GLM and save to analysisData.glm
scanData = getScanData_GLM(thisView,glmInfo.analysisNames_Scans,glmInfo.analysisNames_Groups,glmInfo.groupNames);

%% get data from ROIs
% save so don't need to load again
% Set group outside of script
roiData = script_getROIdata(thisView,scanData.scan_GLMdata,glmInfo.analysisBaseNames_Scans,Info.ROInames,glmInfo.analysisScanNum,'GLM');

%% Convert data to flatmap space

% export scan data
for iScan = 1:glmInfo.nScans
    for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
            analysisName = [glmInfo.analysisBaseNames_Scans{iAnal}, '_Scan_' mat2str(iScan)];
        for iSide = 1:length(subjectInfo.flatmapNames)
            thisView = script_covertData2FlatmapSpace(thisView,'MotionComp',analysisName,iScan,[],subjectInfo.flatmapNames{iSide});
        end
    end
end

% average overdepth
for iSide = 1:length(subjectInfo.flatmapNames)
    thisView = script_averageAcrossDepths(thisView,[],[subjectInfo.flatmapNames{iSide}, 'Volume']);
end

%% export data for concat groups
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},glmInfo.analysisNames_Groups{iAnal},[],[39, 40],subjectInfo.flatmapNames{iSide});
        end
    end
end

%% get data from Overlays
% create names to get data from overlays and save using structure side.Group.anal.data{iScan}
q = char(39);
for iScan = 1:glmInfo.nScans
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
            eval([Info.Sides{iSide}, '.scans.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',conditionNames{iAnal},iScan);'])
        end
    end
end
%% restrict by roi
% script 
% get data from left flatmap (set: left roi, left group, left base) save this data to left struct
% loop - save (side), selet:roi (side) and base (side), get data from: side.scans.anal.overlayData.data
% use eval
%% get data from ROIs
% save so don't need to load again
% Set group outside of script
data_rightROI = script_getROIdata(thisView,Right,glmInfo.analysisBaseNames_Scans,{'FlatRightAC'},glmInfo.analysisScanNum,'overlays');
% change to flat roi names

%% perform ROI analysis
% save so don't need to load again
% compare binning for glm to averaging betas
roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');


roiAnalysis = script_ROIAnalysis(data_rightROI,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays');

%% save data
saveLocation = cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
saveName = 'ROIanalysis.mat';
% navigate to correct directory
save(saveName,'roiAnalysis')
% save data to disk
% load later for group analysis

script_GroupAnalysis(subjects)
% load subject ROI analysis from disk




