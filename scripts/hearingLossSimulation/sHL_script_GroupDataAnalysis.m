function groupData = sHL_script_GroupDataAnalysis
% A script to load individual subject data and perform group level data analysis
%
%   usage: sHL_script_GroupDataAnalysis
%      by: Ben Gurer
%    date: 17/11/2017
% purpose: perform group level data analysis data from  hearing loss simulation study
%   input: n/a


%% get study info
[stimInfo, glmInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects
nSubjects = 3;

%% get subject info
for iSub = 1:nSubjects
subjectInfo(iSub) = get_SubjectInfo_sHL(iSub);
% Subject ID
% flatmap names

%% load data
% move to data location then load and rename or save to subject name struct
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo(iSub).subjectID));
load([subjectInfo.subjectID '_ROIdata.mat']);
end

%% analysis to compare:
% GLM
% ROI average beta weight
% ROI average tuning curves

% pRF
% voxel pCF comparisions (scatter plot)
% voxel pCF distribution

%% pre process data
% normalise 
%   divide by maximium value in condition A (normal hearing)

%% Plot data

end