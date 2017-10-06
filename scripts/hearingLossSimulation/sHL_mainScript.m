sHL_mainScript
% simulate hearing loss
% hearingLossSimulation_mainScript

% setup study direction (once)
sHL_createStudyDirectory

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

studyParamsData = sHL_setupStudyParams;
% stimulus info
% condition names

glmData = script_glmAnalysis(thisView);

thisView = script_importAnatomy(thisView);

pRFdata = script_pRFAnalysis(thisView);

thisView = script_flatMapAnalysis(thisView,flatmapNames);

roiData = script_ROIAnalysis(thisView);
% save data to disk
% load later for group analysis


script_GroupAnalysis(subjects)
% load subject ROI analysis from disk




