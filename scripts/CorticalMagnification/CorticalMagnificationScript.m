%% export files with the same name and different scan number
% make cell array with subject info and then reference that

%% RUN FREESURFER FIRST

iSubj = 1;

epiDims = [128 128 24 361]; % dims of contin
refScan(iSubj) = 15; % scan before t2 structural

% dataDir = '/home/beng/data'
dataDir = '/home/beng/data';
studyDir = 'CorticalMagnification';
subjects{1} = '03644_012';
niftiBaseName{1} = 'pRFpilot2_';
T2star{1} = '16';
wholeheadPSIR = 'xx';
distCorrectionRef = {'17','18'};
scanList = [13,15,19,22];
freeSurferName{1} = 'kkPSIR_reorient_p7';
% distCorrectionRef = {'20','21'};


cd(['~/data/scanner/' subjects{iSubj}])

system('ptoa -f -q -nii *.PAR')
% or try
!ptoa -f -q -nii *.PAR



mkdir(fullfile(dataDir,studyDir,subjects{iSubj}));
cd(fullfile(dataDir,studyDir,subjects{iSubj}));

mkdir('Etc')
mkdir('Distorted')
mkdir('Anatomy')
mkdir('Raw')
mkdir('Raw/TSeries')
% mkdir ('T2_star')

%% Move whole head PSIR
mkdir(fullfile(dataDir,'Anatomy','originals',subjects{iSubj}));
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} wholeheadPSIR '*.nii']),fullfile(dataDir,'Anatomy','originals',subjects{iSubj}));

%% Move in-plane T2 star
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} T2star{iSubj} '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Anatomy'))

%% Move the functional scans
movefile(fullfile(dataDir,'scanner',subjects{iSubj},'*.nii'),fullfile(dataDir,studyDir,subjects{iSubj},'Distorted'))

%% if distortion correction TEs are not interleaved
cd('Distorted')
if length(distCorrectionRef) == 2
    for i = 1:2
        system(['fslsplit ' niftiBaseName{iSubj} distCorrectionRef{i} '_1_modulus ' niftiBaseName{iSubj} distCorrectionRef{i} '_modulus ']);
        system(['fslsplit ' niftiBaseName{iSubj} distCorrectionRef{i} '_1_phase ' niftiBaseName{iSubj} distCorrectionRef{i} '_phase ']);
    end
    fslmergeArgMod = [];
    fslmergeArgPhase = [];
    for i = 0:4
        fslmergeArgMod = [fslmergeArgMod niftiBaseName{iSubj} distCorrectionRef{1} '_modulus' '000' num2str(i) ' ' niftiBaseName{iSubj} distCorrectionRef{2} '_modulus' '000' num2str(i) ' '];
        fslmergeArgPhase = [fslmergeArgPhase niftiBaseName{iSubj} distCorrectionRef{1} '_phase' '000' num2str(i) ' ' niftiBaseName{iSubj} distCorrectionRef{2} '_phase' '000' num2str(i) ' '];
        
    end
    system(['fslmerge -t ' niftiBaseName{iSubj} distCorrectionRef{1} '_' distCorrectionRef{2} '_modulus ' fslmergeArgMod]);
    system(['fslmerge -t ' niftiBaseName{iSubj} distCorrectionRef{1} '_' distCorrectionRef{2} '_phase ' fslmergeArgPhase]);
end

%% take first fram
system(['fslroi ' niftiBaseName{iSubj} distCorrectionRef{1} '_1_modulus ' niftiBaseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame 0 1']);
%%$ Skull strip
system(['bet ' niftiBaseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame ' niftiBaseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame_brain -f .1 -Z'],'-echo');

if check == 1
   system(['fslview '  niftiBaseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame ' niftiBaseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame_brain']);
     %% find wait for keyboard or use input
end

%% now do distortion correction
keyboard
cd('..')
%% copy distortion corrected tseries to raw time series folder
system('cp Distorted/dynB0map/*dynMod_U.nii Raw/TSeries/');


% crop last frame of reference EPI
!mkdir FNIRT
cd Raw/TSeries/


system(['fslroi ' niftiBaseName{iSubj} num2str(refScan(iSubj)) '_1_modulus_dynMod_U.nii ' lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);

!mv lastFrameEPI.nii ../../FNIRT
cd ../../

% crop T2*

cd Anatomy/
T2starFile = [niftiBaseName{iSubj} T2star{iSubj} '_1_modulus'];
system(sprintf('fslroi %s %s_crop 18 348 18 348 2 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));
system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_stripped -f .2 -Z'])

cd ../

%run FNIRT
cd FNIRT
system(sprintf('fnirtEpi2T2star lastFrameEPI %s_crop -separateFLIRT',T2starFile));

mrAlign
keyboard
% manually align the go to
% Compute Alignment->Advanced Alignment Menu
% then select Reverse Contrast (T2*)
% Click 'Compute Coarse Aligment' then
% Click 'Compute fine Aligment'

% once done
% save alignment to file and files (structural T2* images of the SAME size)

% Reload source as destination
% load lastFrameEPI as source
% Click Manual Alignment-> Set Alignment to Identity
% Save alignment to files (all functional data, including distorted)
%

%% align functional data to high resolution t2* inplane structural image

% for i = 1:length(scanList)
% alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],fullfile(dataDir,studyDir,subjects{iSubj},'Raw/TSeries',[niftiBaseName{iSubj} num2str(scanList(i)) '_1_modulus_dynMod_U.nii']));
% end
lastFrameNames = struct2cell(dir('lastFrameEPI*.nii'));
lastFrameNames = lastFrameNames(1,:);
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],lastFrameNames);

functionalNames = struct2cell(dir('../Raw/TSeries/*.nii'));
functionalNames = functionalNames(1,:);
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],functionalNames);

cd ..


%Pre-process functional

[sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1'); % looks in Raw/Tseries to find out how many scans there are
sessionParams.subject = subjects{iSubj};
sessionParams.description = 'CorticalMagnification';
sessionParams.operator = 'bg';

%% This is just for this pilot - change when we know what we are doing
groupParams.description([1,3]) = {'Sparse, Run 1','Sparse, Run 2'};
groupParams.description([2,4]) = {'Continuous, Run 1','Continuous, Run 2'};
nScans = length(groupParams.name);
% for iScan = 1:nScans
%   groupParams.description{iScan} = ['Run ' num2str(iScan)];
% end

mrInit(sessionParams,groupParams,'makeReadme=0');

%motion correction
thisView = newView;
refScanNum = viewGet(thisView,'scannum',sprintf('%s%02d_1_modulus_dynMod_U.nii',niftiBaseName{iSubj},refScan(iSubj)));
[thisView, motionCompParams] = motionComp(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(1:nScans)]);
motionCompParams.baseFrame='last';
motionCompParams.baseScan = refScanNum;
[thisView, motionCompParams] = motionComp(thisView,motionCompParams);

%concatenation
thisView = viewSet(thisView,'curGroup','MotionComp');
[thisView, concatParams] = concatTSeries(thisView,[],'defaultParams=1',['scanList=' mat2str(1:nScans)]);


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

%load reference EPI as anatomy
thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT'));
%load skull-stripped  EPI as overlay
[thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(niftiBaseName, subjects{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

mrSaveView(thisView);
deleteView(thisView);

save('preProcessParams.mat','motionCompParams','concatParams');



mrLoadRet


% THE REST REQUIRES THE FLAT MAPS
keyboard

%IMPORT  FREESURFER SURFACE 
cd(fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj}));
mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[336 336 336]'); % this will depend on the resolution of the PSIR

%apply FNIRT warping coefficient to surfaces
cd(fullfile(dataDir,studyDir,subjects{iSubj}));
fslApplyWarpSurfOFF(fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/',[niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_resampled_2_lastFrameEPI_warpcoef.nii']),...
                    fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/','lastFrameEPI.nii'),...
                    fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj},'surfRelax',[freeSurferName{iSubj} '_mprage_pp.nii']),...
                    subjects{iSubj});
                
% add how to do freesurfer stuff
