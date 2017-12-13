function [thisView, concatedata] = script_setupmrLoadRet(Info,subjectInfo,glmInfo)
% function to script creating a new mrLoadRet view
% creates view
% performs motion correction on scans
% concatinates scans into groups


%% Set up mrTools mrLoadRet
[sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1'); % looks in Raw/Tseries to find out how many scans there are
sessionParams.subject = subjectInfo.subjectID;
sessionParams.description = Info.studyDir;
sessionParams.operator = 'bg';

% define in subject info
groupParams.description(subjectInfo.conditionOrder{2}) = {'Hearing Loss Simulation, Run 1','Hearing Loss Simulation, Run 2'};
groupParams.description(subjectInfo.conditionOrder{1}) = {'Normal Hearing, Run 1','Normal Hearing, Run 2'};
nScans = length(groupParams.name);

mrInit(sessionParams,groupParams,'makeReadme=0');

% Motion correction
thisView = newView;
refScanFileName = [subjectInfo.niftiBaseName 'WIP_73DYN_fMRI_02_' num2str(subjectInfo.refScan) '*.nii'];   
refScanNum = viewGet(thisView,'scannum',refScanFileName);
[thisView, motionCompParams] = motionComp(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(1:nScans)]);
motionCompParams.baseFrame='last';
motionCompParams.baseScan = refScanNum;
[thisView, motionCompParams] = motionComp(thisView,motionCompParams);

% Concatenation of Normal Hearing data
thisView = viewSet(thisView,'curGroup','MotionComp');
params_ConcatenationNH = getConcatParams_withNewGroupName(thisView,'ConcatenationNH','defaultParams=1',['scanList=' mat2str(subjectInfo.conditionOrder{1})]);
[thisView, params_ConcatenationNH] = concatTSeries(thisView,params_ConcatenationNH);

% Concatenation of Hearing Loss Simulation data
thisView = viewSet(thisView,'curGroup','MotionComp');
params_ConcatenationHLsim = getConcatParams_withNewGroupName(thisView,'ConcatenationHLsim','defaultParams=1',['scanList=' mat2str(subjectInfo.conditionOrder{2})]);
[thisView, params_ConcatenationHLsim] = concatTSeries(thisView,params_ConcatenationHLsim);

% link stim files to scans
system(sprintf('cp %s/*.txt Etc/',fullfile(Info.dataDir,'scanner',subjectInfo.subjectID,'logFiles')));
cd Etc/
logFiles = dir('*.txt');
logToMylogAdaptation([],{logFiles(:).name});
cd ..
logFiles = dir('Etc/*.mylog.mat');
logFiles = {logFiles(:).name};
for iFile = 1:length(logFiles)
    fprintf(1,['Linking ' logFiles{iFile} ' to Group Raw, scan ' num2str(viewGet(thisView,'tseriesFile',iFile,1)) '\n']);
    thisView = viewSet(thisView,'stimfilename',logFiles{iFile}, iFile,1);
end

save('preProcessParams.mat','motionCompParams','params_ConcatenationNH','params_ConcatenationHLsim');

% save view and quit
mrSaveView(thisView);

end