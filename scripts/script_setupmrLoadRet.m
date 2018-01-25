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
params_ConcatenationNH = getConcatParams_withNewGroupName(thisView,glmInfo.groupNames{1},'defaultParams=1',['scanList=' mat2str(subjectInfo.conditionOrder{1})]);
[thisView, params_ConcatenationNH] = concatTSeries(thisView,params_ConcatenationNH);

% Concatenation of Hearing Loss Simulation data
thisView = viewSet(thisView,'curGroup','MotionComp');
params_ConcatenationHLsim = getConcatParams_withNewGroupName(thisView,glmInfo.groupNames{2},'defaultParams=1',['scanList=' mat2str(subjectInfo.conditionOrder{2})]);
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


%% smooth and create new group
if Info.smooth == 1
    % cd to motioncomp group
    % get file names
    % loop
    % run:
    % fslmaths original.nii -kernel gauss 2.1233226 -fmean smoothed.nii
    % create new group
    % move scans to new folder
    % add scans
    % concatenate runs
    
    % move file directory to MotionComp group
    cd MotionComp/TSeries/
    
    % change mrTools view to MotionComp group
%     thisView = viewSet(thisView,'curGroup','MotionComp');
    
    % get all file names from directory
    scanFiles = [];
    scanFiles = dir('*.nii');
    scanFiles = {scanFiles(:).name};
    
    % Smooth scans using fslmaths
    sigma = Info.smoothingFWHM./ 2.355;
%     sigma = 2;
    for iFile = 1:length(scanFiles)
        % fslmaths original.nii -kernel gauss 2.1233226 -fmean smoothed.nii
%         system(['fslmaths ' scanFiles{iFile} ' -kernel gauss ' num2str(sigma) ' -fmean' [ '_smoothed_' scanFiles{iFile} ] '.nii']);
    system(['fslmaths ' scanFiles{iFile} ' -s ' num2str(sigma) ' ' [ 'smoothed_' scanFiles{iFile} ] '.nii']);
  
    end
    smoothedScanFiles = [];
%     smoothedScanFiles = dir('-fmean*.nii');
    
    smoothedScanFiles = dir('smoothed_*.nii');
    smoothedScanFiles = {smoothedScanFiles(:).name};
    % move back to subject root directory
    cd ../..
    
    thisView = getMLRView;
    % create 'Smoothed' group in mrTools view
%     groupName = ['Smoothed_fwhm' num2str(round(gaus_fwhm))];
    groupName = Info.smoothingName;
    thisView = viewSet(thisView,'newGroup',groupName);
        
    % change mrTools view to 'Smoothed' group
    thisView = viewSet(thisView,'curGroup',groupName);
    
    
   smoothedGroupNum = viewGet(thisView,'groupnum',groupName);
    
    % import scans and link log files
    
    logFiles = [];
    logFiles = dir('Etc/*.mylog.mat');
    logFiles = {logFiles(:).name};
    % change to importTSeires from motion comp and then delete ffrom
    % motioncomp after
    for iFile = 1:length(logFiles)
        thisView = importTSeries(thisView,[],'defaultParams=1',['pathname=' fullfile('MotionComp','TSeries',smoothedScanFiles{iFile})]);
        fprintf(1,['Linking ' logFiles{iFile} ' to Group ' groupName ', scan ' num2str(viewGet(thisView,'tseriesFile',iFile,1)) '\n']);
        thisView = viewSet(thisView,'stimfilename',logFiles{iFile}, iFile,smoothedGroupNum);
    end
 
    % move smoothed data to Smoothed group directory
    !rm -f MotionComp/TSeries/smoothed_*.nii
    
    for iGroup = 1:length(glmInfo.groupNames)
        params_Concatenation = getConcatParams_withNewGroupName(thisView,glmInfo.groupNames{iGroup},'defaultParams=1',['scanList=' mat2str(subjectInfo.conditionOrder{iGroup})]);
        [thisView, params_Concatenation] = concatTSeries(thisView,params_Concatenation);
    end
    
end

% save view and quit
mrSaveView(thisView);

refreshMLRDisplay(thisView.viewNum);

end