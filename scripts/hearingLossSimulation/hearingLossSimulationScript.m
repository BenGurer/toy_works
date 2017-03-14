iSubj = 1;

epiDims = [128 128 24 73]; % dims of functional scans

% add check from startup

if ispc
    dataDir = 'N:/data';
elseif isunix
    dataDir = '/home/beng/data';
end
studyDir = 'hearingLossSimulation';

sides = {'left','right'};
Sides = {'Left','Right'};

% Subject info
subjects{1} = '02344_034';
niftiBaseName{1} = 'Sparse_02344_SPARSE_64dYN_1.5mmiso_TE40_SENSE_';
wholeheadMPRAGE{1} = 'Sparse_02344_MPRAGE_SENSE_13_1';
freeSurferName{1} = '02344_034';
T2star{1} = 'Sparse_02344_High_res_t2__SENSE_12_1';
refScan{1} = '11'; % scan before t2 structural

% distCorrectionRefSparse{1} = {'17','18'};
% distCorrectionRefCont{2} = {'20','21'};
% % scanList = [13,15,19,22];
% distCorrectionRef = {'20','21'};

cd([dataDir '/scanner/' subjects{iSubj}]) 
 
!ptoa -f -q -nii *.PAR

if ~isempty(wholeheadMPRAGE{iSubj})
    % Move whole head PSIR
    mkdir(fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    movefile(fullfile(dataDir,'scanner',subjects{iSubj},[wholeheadMPRAGE{iSubj} '*.nii']),fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    cd(fullfile(dataDir,'Anatomy/originals/',freeSurferName{iSubj}));        
    %RUN recon-all DIRECTLY IN TERMINAL 
    fprintf(['recon-all -subjid ' freeSurferName{iSubj} ' -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' wholeheadMPRAGE{iSubj} '.nii' ' -all']);
end

mkdir(fullfile(dataDir,studyDir,subjects{iSubj}));
cd(fullfile(dataDir,studyDir,subjects{iSubj}));

mkdir('Etc')
mkdir('Anatomy')
mkdir('Raw')
mkdir('Raw/TSeries')
!mkdir FNIRT

%% Move in-plane T2 star
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[T2star{iSubj} '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Anatomy'))
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Raw/TSeries'))

%% check for filenames below 10 and add a zero before the number - make sure linking of stim files works later
cd Raw/TSeries/
scanFiles = dir;
for id = 1:length(scanFiles)
    str = scanFiles(id).name;
    strParts = strsplit(str,'_');
    if length(strParts) > 1
    checkScanNum = char(strParts(8));
%         numStrParts = length(strParts);
        if 10>str2num(checkScanNum) && numel(checkScanNum)==1
            newName = [subjects{iSubj} '_0' char(strParts(8)) '.nii']
            movefile(scanFiles(id).name,newName);
        else
                        newName = [subjects{iSubj} '_' char(strParts(8)) '.nii']
            movefile(scanFiles(id).name,newName);
        end
   
    end
end
% crop last frame of reference EPI

system(['fslroi ' subjects{iSubj} '_' refScan{iSubj} '.nii lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);

!mv lastFrameEPI.nii ../../FNIRT
cd ../../

% crop T2*
cd Anatomy/
mrAlign
keyboard
% align hi res t2* to MPRAGE
% Load MPRAGE as destination and ' Set base cooridinate frame'
% Load hi res t2* as source
% manually align then go to
% Compute Alignment->Advanced Alignment Menu
% then select Reverse Contrast (T2*)
% Click 'Compute Coarse Aligment' then
% Click 'Compute fine Aligment'
% save alignment
% doing it here means cropped versions with have the correct alignement -
% skull stripping removes the s-form matrix so using set alignment to
% identity and save (t2*_stripped to t2*)

T2starFile = [T2star{iSubj} '_modulus'];
T2StarDims_Index = [400 400 29];
EPIDims_Index = [128 128 24];
T2StarDims_mm = [0.55 0.55 1.5];
EPIDims_mm = [1.359375  1.359375 1.5];
% 400 400 29
% 0.55 0.55 1.5
% 220 220 1
% 128 128 24
% 1.359375  1.359375 1.5
% 174 mm diff = 46 mm
% 316
% 46/0.55 = 83 /2 = 42
% fslroi input = (first index) (number of frames)
% system(sprintf('fslroi %s %s_crop 22 316 22 316 3 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('fslroi %s %s_crop 42 316 42 316 3 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data

% system(sprintf('fslroi %s %s_crop 18 348 18 348 2 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));
% system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_stripped -f .2 -Z'])
system(['bet ' T2star{iSubj} '_modulus ' T2star{iSubj} '_modulus_stripped -f .2 -Z'])

mrAlign
keyboard
% skull stripping removes the s-form matrix so using set alignment to
% identity and save (t2*_stripped to t2*)

cd ../FNIRT

%run FNIRT
system(sprintf('fnirtEpi2T2star lastFrameEPI %s_crop -separateFLIRT',T2starFile));

lastFrameNames = struct2cell(dir('lastFrameEPI*.nii'));
lastFrameNames = lastFrameNames(1,:);
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],lastFrameNames);

functionalNames = struct2cell(dir('../Raw/TSeries/*.nii'));
functionalNames = functionalNames(1,:);
for i = 1:length(functionalNames)
    functionalNames{i} = ['../Raw/TSeries/' functionalNames{i}];
end
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],functionalNames);

cd ../


%Pre-process functional#


%% Set up mrTools mrLoadRet

[sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1'); % looks in Raw/Tseries to find out how many scans there are
sessionParams.subject = subjects{iSubj};
sessionParams.description = studyDir;
sessionParams.operator = 'bg';

groupParams.description([1,3]) = {'Hearing Loss Simulation, Run 1','Hearing Loss Simulation, Run 2'};
groupParams.description([2,4]) = {'Normal Hearing, Run 1','Normal Hearing, Run 2'};
nScans = length(groupParams.name);

mrInit(sessionParams,groupParams,'makeReadme=0');

%motion correction
thisView = newView;
refScanNum = viewGet(thisView,'scannum',sprintf('%s_%s.nii',subjects{iSubj},refScan{iSubj}));
[thisView, motionCompParams] = motionComp(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(1:nScans)]);
motionCompParams.baseFrame='last';
motionCompParams.baseScan = refScanNum;
[thisView, motionCompParams] = motionComp(thisView,motionCompParams);

%concatenation of Hearing Loss Simulation data
thisView = viewSet(thisView,'curGroup','MotionComp');
params = getConcatParams_withNewGroupName(thisView,'ConcatenationHLsim','defaultParams=1',['scanList=' mat2str([1 3])]);
[thisView, concatParamsSparse] = concatTSeries(thisView,params);

%concatenation of Normal Hearing data
thisView = viewSet(thisView,'curGroup','MotionComp');
params = getConcatParams_withNewGroupName(thisView,'ConcatenationNH','defaultParams=1',['scanList=' mat2str([2 4])]);
[thisView, concatParamsCont] = concatTSeries(thisView,params);

%GLM analysis
system(sprintf('cp %s/*.txt Etc/',fullfile(dataDir,'scanner',subjects{iSubj},'logFiles')));
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

thisView = viewSet(thisView,'curGroup','ConcatenationHLsim');
[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.hrfModel = 'hrfBoxcar';
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_BoxCar';
glmParams.hrfParams.description = 'GLM Box Car -HLsim Concat';
glmParams.hrfParams.delayS =  2.5;
glmParams.hrfParams.durationS = 2.5;
glmParams.scanParams{1}.stimDurationMode = 'fromFile';
glmParams.scanParams{1}.supersamplingMode =  'Set value';
glmParams.scanParams{1}.designSupersampling = 3;
glmParams.scanParams{1}.acquisitionDelay = .75;
glmParams.numberFtests = 1;
glmParams.fTestNames{1} = 'fTest 1';
glmParams.numberContrasts = 0;
glmParams.parametricTests = 0;
glmParams.outputEstimatesAsOverlays = 1; 
[thisView, glmParams] = glmAnalysis(thisView,glmParams);

%Tonotopy analysis
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);

params.combineFunction='weightedMeanStd';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);

%% save analysis
saveAnalysis(thisView,'GLM_BoxCar')

thisView = viewSet(thisView,'curGroup','ConcatenationNH');
[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.hrfModel = 'hrfBoxcar';
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_BoxCar';
glmParams.hrfParams.description = 'GLM Box Car -NH Concat';
glmParams.hrfParams.delayS =  2.5;
glmParams.hrfParams.durationS = 2.5;
glmParams.scanParams{1}.stimDurationMode = 'fromFile';
glmParams.scanParams{1}.supersamplingMode =  'Set value';
glmParams.scanParams{1}.designSupersampling = 3;
glmParams.scanParams{1}.acquisitionDelay = .75;
glmParams.numberFtests = 1;
glmParams.fTestNames{1} = 'fTest 1';
glmParams.numberContrasts = 0;
glmParams.parametricTests = 0;
glmParams.outputEstimatesAsOverlays = 1; 
[thisView, glmParams] = glmAnalysis(thisView,glmParams);

%Tonotopy analysis
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);

params.combineFunction='weightedMeanStd';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);


%% save analysis
saveAnalysis(thisView,'GLM_BoxCar')

% %load reference EPI as anatomy
% thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT'));
%load skull-stripped  EPI as overlay
% [thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

% save('preProcessParams.mat','motionCompParams','concatParams');
save('preProcessParams.mat','motionCompParams');

% script running glm on each scan individually and performing split
% analysis


%IMPORT  FREESURFER SURFACE
cd(fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj}));
mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[256 256 162]'); % 256 256 162

%apply FNIRT warping coefficient to surfaces
cd(fullfile(dataDir,studyDir,subjects{iSubj}));
fslApplyWarpSurfOFF(fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/',[niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_resampled_2_lastFrameEPI_warpcoef.nii']),...
    fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/','lastFrameEPI.nii'),...
    fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj},'surfRelax',[freeSurferName{iSubj} '_mprage_pp.nii']),...
    subjects{iSubj});

%import surfaces (this step has not been made scriptable yet) 
for iSide=1:2
  base = importSurfaceOFF(fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj},'surfRelax',...
   [freeSurferName{iSubj} '_' sides{iSide} '_GM.off']));
  thisView = viewSet(thisView, 'newbase', base);
  thisView = viewSet(thisView,'corticalDepth',[0.2 0.8]);
end


% save view and quit
mrSaveView(thisView);
deleteView(thisView);


mrLoadRet
% wait until it loads
thisView = getMLRView;

flatmapInfo{1} = {'80_131_81_Rad60', '181_127_85_Rad60'};
concatenationGroup = {'ConcatenationHLsim', 'ConcatenationNH'};
functionalAnalysis = {'GLM_BoxCar'};
mainOverlays(1) = 38;
ROInames = {'LeftPAC','RightPAC'};

% analysisType = viewGet(thisView,'analysisType');
roiNum = 2;
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
% need to export to new group because over rights / need one for each
% flatmap/hemisphere
for iGroup = 1:length(concatenationGroup)
for iAnalysis = 1:length(functionalAnalysis)
for iSide=1:2
  % gradient reversals
  thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
  thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis{iAnalysis}));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['x' freeSurferName{iSubj} '_' sides{iSide} '_WM_Flat_' flatmapInfo{iSubj}{iSide}]));
  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(mainOverlays(iSubj,1))]);
  params.combineFunction='gradientReversal';
  params.additionalArgs = '[18 18 21]';
  params.baseSpaceInterp = 'linear';
  params.nOutputOverlays=7;
  params.baseSpace = 1;
%   params.exportToNewGroup=1;
  [thisView,params] = combineTransformOverlays(thisView,params);
  curOverlay=viewGet(thisView,'curOverlay');
  thisView = viewSet(thisView,'overlayMin',15,curOverlay-1);
  thisView = viewSet(thisView,'overlayMax',180,curOverlay-1);
  thisView = viewSet(thisView,'overlaycolorRange',[45 180],curOverlay-1);
  thisView = viewSet(thisView,'overlayMax',75);
  thisView = viewSet(thisView,'overlaycolorRange',[0 90]);
end
end
end

for iGroup = 1:length(concatenationGroup)
for iAnalysis = 1:length(functionalAnalysis)
    %roi loop here
  thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
  thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis{iAnalysis}));
  
  analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',functionalAnalysis{iAnalysis}));
  glmData{iGroup} = analysisData.d{1};
  analysisParams{iGroup} = analysisData.params;
  r2data = analysisData.overlays(1).data{1};
  roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
  %get ROI estimates
  volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
  [e_contGroups{iGroup},volumeIndices] = getEstimates(glmData{iGroup} ,analysisParams{iGroup} ,volumeIndices');
  nVoxels = length(volumeIndices);
end
end
[ROIbetasSum_groups{1}, ROIStesSum_groups{1}] = plotROIav_GLMBetaEstimates_SplitRuns(e_contGroups{1},e_contGroups{2});
betasCondA = squeeze(e_contGroups{1}.betas);
betasCondB = squeeze(e_contGroups{2}.betas);

figure
plot(sum(betasCondA,2))
figure
plot(sum(betasCondB,2))

stimBins = 8;
for iScan = 1:nScans
thisView = viewSet(thisView,'curGroup','MotionComp',['curScan=' mat2str(iScan)]);
[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
glmParams.hrfModel = 'hrfBoxcar';
glmParams.hrfParams.delayS =  2.5;
glmParams.hrfParams.durationS = 2.5;
% [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
glmParams.hrfParams.description = ['GLM Box Car - Scan ' mat2str(iScan)];
glmParams.scanParams{1, iScan}.preprocess  = 'binStimFreq';
glmParams.numberContrasts = stimBins;
glmParams.numberEVs = stimBins;
glmParams.EVnames = {'1','2','3','4','5','6','7','8'};
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
glmParams.scanParams{iScan}.stimDurationMode = 'fromFile';
glmParams.scanParams{iScan}.supersamplingMode =  'Set value';
glmParams.scanParams{iScan}.designSupersampling = 3;
glmParams.scanParams{iScan}.acquisitionDelay = .75;

glmParams.numberFtests = 1;
glmParams.fTestNames{1} = 'fTest 1';
glmParams.numberContrasts = 0;
glmParams.parametricTests = 0;
glmParams.outputEstimatesAsOverlays = 1; 
glmParams.saveName = ['GLM Box Car - Scan ' mat2str(iScan)];
[thisView, glmParams] = glmAnalysis(thisView,glmParams,['scanList=' mat2str(iScan)]);

%Tonotopy analysis
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:stimBins+1])],['scanList=' mat2str(iScan)]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 stimBins],curOverlay-1);
% 
params.combineFunction='weightedMeanStd';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 stimBins],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 stimBins],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 stimBins*1.25],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 stimBins*1.25],curOverlay);

% save analysis
saveAnalysis(thisView,['GLM_BoxCar - Scan ' mat2str(iScan)])
end


thisView = getMLRView;
ROInames = {'LeftPAC','RightPAC'};

analysisType = viewGet(thisView,'analysisType');
% analysisParams = convertOldGlmParams(viewGet(thisView,'analysisParams'));
roiNum = 2;
% 
%     roiList = selectInList(thisView,'rois');
%     
%   roi = viewGet(thisView,'roi',iRoi);
%   
%   roi{roiNum} = viewGet(thisView,'roi','RightPAC');
  
  roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
  e = cell(1,4);
  thisView = viewSet(thisView,'curgroup','MotionComp');
for iScan = 1:nScans
% for iSide=1:2
%   thisView = viewSet(thisView,'curgroup','MotionComp',['curScan=' mat2str(iScan)]);
  
  thisView = viewSet(thisView,'curScan' ,iScan);
%    refreshMLRDisplay(thisView.viewNum);
  thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',['GLM Box Car - Scan ' mat2str(iScan)]));
%   refreshMLRDisplay(thisView.viewNum);
  
% r2data = viewGet(thisView,'overlaydata',1);

 
analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',['GLM Box Car - Scan ' mat2str(iScan)]));
glmData{iScan} = analysisData.d{iScan};
analysisParams{iScan} = analysisData.params;
r2data = analysisData.overlays(1).data{iScan}; 
% glmData = viewGet(thisView,'d',viewGet(thisView,'analysisNum',['GLM Box Car - Scan ' mat2str(iScan)]))
% glmData = viewGet(thisView,'d',viewGet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',['GLM Box Car - Scan ' mat2str(iScan)])));
      % get roi scan coords
    roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
        %get ROI estimates 
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
%     roiIndices = (r2data(volumeIndices)>r2clip(1)) & (r2data(volumeIndices)<r2clip(2));% & (~isnan(volumeBetas(volumeIndices,1,1)))';
%     volumeIndices = volumeIndices(roiIndices);
    [e{iScan},volumeIndices] = getEstimates(glmData{iScan} ,analysisParams{iScan} ,volumeIndices');
    nVoxels = length(volumeIndices);
%     nTotalVoxels = length(roiIndices);
  
% end
end
[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns(e{1},e{3});

[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns(e{2},e{4});

