%% export files with the same name and different scan number
% make cell array with subject info and then reference that

%% RUN FREESURFER FIRST

iSubj = 1;

epiDims = [128 128 24 361]; % dims of contin
refScan(iSubj) = 15; % scan before t2 structural


mountPoint = '/home/beng/data';
studyName = 'CorticalMagnification';
subjectID{1} = '03644_012';
baseName{1} = 'pRFpilot2_';
inplaneT2star{1} = '16';
wholeheadPSIR = 'xx';
distCorrectionRef = {'17','18'};
% distCorrectionRef = {'20','21'};


cd(['~/data/scanner/' subjectID{iSubj}])

system('ptoa -f -q -nii *.PAR')
% or try
!ptoa -f -q -nii *.PAR



mkdir(fullfile(mountPoint,studyName,subjectID{iSubj}));
cd(fullfile(mountPoint,studyName,subjectID{iSubj}));

mkdir('Etc')
mkdir('Distorted')
mkdir('Anatomy')
mkdir('Raw')
mkdir('Raw/TSeries')
% mkdir ('T2_star')

%% Move whole head PSIR
mkdir(fullfile(mountPoint,'Anatomy','originals',subjectID{iSubj}));
movefile(fullfile(mountPoint,'scanner',subjectID{iSubj},[baseName{iSubj} wholeheadPSIR '*.nii']),fullfile(mountPoint,'Anatomy','originals',subjectID{iSubj}));

%% Move in-plane T2 star
movefile(fullfile(mountPoint,'scanner',subjectID{iSubj},[baseName{iSubj} inplaneT2star{iSubj} '*.nii']),fullfile(mountPoint,studyName,subjectID{iSubj},'Anatomy'))

%% Move the functional scans
movefile(fullfile(mountPoint,'scanner',subjectID{iSubj},'*.nii'),fullfile(mountPoint,studyName,subjectID{iSubj},'Distorted'))

%% if distortion correction TEs are not interleaved
cd('Distorted')
if length(distCorrectionRef) == 2
    for i = 1:2
        system(['fslsplit ' baseName{iSubj} distCorrectionRef{i} '_1_modulus ' baseName{iSubj} distCorrectionRef{i} '_modulus ']);
        system(['fslsplit ' baseName{iSubj} distCorrectionRef{i} '_1_phase ' baseName{iSubj} distCorrectionRef{i} '_phase ']);
    end
    fslmergeArgMod = [];
    fslmergeArgPhase = [];
    for i = 0:4
        fslmergeArgMod = [fslmergeArgMod baseName{iSubj} distCorrectionRef{1} '_modulus' '000' num2str(i) ' ' baseName{iSubj} distCorrectionRef{2} '_modulus' '000' num2str(i) ' '];
        fslmergeArgPhase = [fslmergeArgPhase baseName{iSubj} distCorrectionRef{1} '_phase' '000' num2str(i) ' ' baseName{iSubj} distCorrectionRef{2} '_phase' '000' num2str(i) ' '];
        
    end
    system(['fslmerge -t ' baseName{iSubj} distCorrectionRef{1} '_' distCorrectionRef{2} '_modulus ' fslmergeArgMod]);
    system(['fslmerge -t ' baseName{iSubj} distCorrectionRef{1} '_' distCorrectionRef{2} '_phase ' fslmergeArgPhase]);
end

%% take first fram
system(['fslroi ' baseName{iSubj} distCorrectionRef{1} '_1_modulus ' baseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame 0 1']);
%%$ Skull strip
system(['bet ' baseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame ' baseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame_brain -f .1 -Z'],'-echo');

if check == 1
   system(['fslview '  baseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame ' baseName{iSubj} distCorrectionRef{1} '_modulus_firstFrame_brain']);
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


system(['fslroi ' baseName{iSubj} num2str(refScan(iSubj)) '_1_modulus_dynMod_U.nii ' lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);

!mv lastFrameEPI.nii ../../FNIRT
cd ../../

% crop T2*

cd Anatomy/
T2starFile = [baseName{iSubj} inplaneT2star{iSubj} '_1_modulus'];
system(sprintf('fslroi %s %s_crop 18 348 18 348 2 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));
system(['bet ' baseName{iSubj} inplaneT2star{iSubj} '_1_modulus_crop ' baseName{iSubj} inplaneT2star{iSubj} '_1_modulus_crop_stripped -f .2 -Z'])

cd ../

%run FNIRT
cd FNIRT
system(sprintf('fnirtEpi2T2star lastFrameEPI %s_crop',T2starFile));

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
cd ..


%Pre-process functional

[sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1'); % looks in Raw/Tseries to find out how many scans there are
sessionParams.subject = subjectID{iSubj};
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
refScanNum = viewGet(thisView,'scannum',sprintf('%s%02d_1_modulus_dynMod_U.nii',baseName{iSubj},refScan(iSubj)));
[thisView, motionCompParams] = motionComp(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(1:nScans)]);
motionCompParams.baseFrame='last';
motionCompParams.baseScan = refScanNum;
[thisView, motionCompParams] = motionComp(thisView,motionCompParams);

%concatenation
thisView = viewSet(thisView,'curGroup','MotionComp');
% [thisView, concatParams] = concatTSeries(thisView,[],'defaultParams=1',['scanList=' mat2str(1) mat2str(3)]);
% [thisView, concatParams] = concatTSeries(thisView,[],'defaultParams=1',['scanList=' mat2str(2) mat2str(4)]);

%load reference EPI as anatomy
thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(baseName, subjectID{iSubj},'FNIRT'));
%load skull-stripped  EPI as overlay
[thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(baseName, subjectID{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

mrSaveView(thisView);
deleteView(thisView);

save('preProcessParams.mat','motionCompParams');

mrLoadRet
%%%%%%%%%%%%%% DEFINE AUDITORY CORTEX ROIs BASED ON F-TEST AND BRAIN MASK
keyboard

