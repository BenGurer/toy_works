function CM_groupAnalysisScript

% get study info
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;

% use is pc to set data directory - could do in cm_setupStduyparams
% Info.dataDir = '/Volumes/data_PSY/data';
Info.dataDir = 'E:\data';
q = char(39);

%% Load subject data
groupData = struct;
% for iSub = 1:8
for iSub = 5
% Get subject info
subjectInfo = get_SubjectInfo_CM(iSub);
% Subject ID, flatmap names
saveName = [subjectInfo.subjectID '_data.mat'];
% move to subject folder
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
% load data
groupData(iSub) = load(saveName);
end

%% Tidy data
% Transform data into tidyVerse
% Row: observation = voxel
% Column: Variable = everything else
% side, group(acquisiton, run), analysis, estimation method, estimate

%% get GLM data
% group

% scans


%% get pRF data

%% Estimating HRF
% average HRF estiamte


%% Compare Analysis



%% Compare Acquisiton
% Sparses Vs Continuous
% use best analysis
% Calculate
% Visualize






end
