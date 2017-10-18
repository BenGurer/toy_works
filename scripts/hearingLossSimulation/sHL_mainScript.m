%% sHL_mainScript - scripted analysis for Simulated Hearing Loss (sHL) study 2017
% simulate hearing loss

%% run this function when returning for analysis?
iSub = 2;

thisView = getMLRView;

% setup study direction (once)
sHL_createStudyDirectory

[stimInfo, glmInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

% define and get subject info
% either cell for each subject or return single subject info - prob the later
subjectInfo = get_SubjectInfo_sHL(iSub);
% Subject ID
% flatmap names


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

thisView = script_importAnatomy(thisView);
% load in:
% reference EPI
% High resolution in-plane T2*
% surfaces
% create flatmaps

[thisView, pRFdata] = script_pRFAnalysis(thisView);

thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);

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
roiData = script_getROIdata(thisView,scanData.scan_GLMdata,glmInfo.analysisBaseNames_Scans,Info.ROInames,glmInfo.analysisScanNum);

%% perform ROI analysis
% save so don't need to load again
roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum);

%% save data
saveName = [subjectID, 'ROIanalysis'];
% navigate to correct dirtory
save(saveName,roiAnalysis)
% save data to disk
% load later for group analysis


script_GroupAnalysis(subjects)
% load subject ROI analysis from disk




