sHL_mainScript
% simulate hearing loss
% hearingLossSimulation_mainScript

% setup study direction (once)
sHL_createStudyDirectory

studyParamsData = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

% define and get subject info
sHL_getSubjectInfo(iSub)
% Subject ID
% flatmap names

sHL_organiseData(studyDirectory,subject)
% per subject

preprocessData = sHL_preprocess;
% distortion correct
% linear alignment
% non-linear alignment

[thisView, concatedate] = script_setupmrLoadRet(thisView,concatinationNames);
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

thisView = script_flatMapAnalysis(thisView,flatmapNames);

%% Create subject data storage structure

%% get data from scans and groups
% save so don't need to load again

%% get data from ROIs
% save so don't need to load again

%% perform ROI analysis
% save so don't need to load again

%% save data

%% new section or script for group analysis

roiData = script_ROIAnalysis(thisView);
% save data to disk
% load later for group analysis


script_GroupAnalysis(subjects)
% load subject ROI analysis from disk




