%% export files with the same name and different scan number
% make cell array with subject info and then reference that

%% RUN FREESURFER FIRST

iSubj = 5;

epiDims = [128 128 24 361]; % dims of contin

% add check from startup

if ispc
    dataDir = 'N:/data';
elseif isunix
    dataDir = '/home/beng/data';
end
studyDir = 'CorticalMagnification';


% Subject info
subjects{1} = '03644_012';
niftiBaseName{1} = 'pRFpilot2_';
T2star{1} = '16';
refScan(1) = 15; % scan before t2 structural
wholeheadPSIR{1} = [];
distCorrectionRefSparse{1} = {'17','18'};
distCorrectionRefCont{2} = {'20','21'};
% scanList = [13,15,19,22];
freeSurferName{1} = 'kkPSIR_reorient_p7';
% distCorrectionRef = {'20','21'};

subjects{2} = '12013_001';
niftiBaseName{2} = 'cm_12013_001_';
psirNiftiBaseName{2} = 'cm_12013_001';
T2star{2} = '13';
refScan(2) = 8; % scan before t2 structural
wholeheadPSIR{2} = '16';
distCorrectionRefSparse{2} = {'9','10'};
distCorrectionRefCont{2} = {'11','12'};
freeSurferName{2} = '12013_001';
sparseScans{2} =  {'7','14'};
contScans{2} =  {'8','15'};

subjects{3} = '12022_001';
niftiBaseName{3} = 'cm_12022_001_';
psirNiftiBaseName{3} = 'cm_12022_001';
T2star{3} = '13';
refScan(3) = 8; % scan before t2 structural
wholeheadPSIR{3} = '17';
distCorrectionRefSparse{3} = {'9','10'};
distCorrectionRefCont{3} = {'11','12'};
freeSurferName{3} = '12022_001';
sparseScans{3} =  {'7','15'};
contScans{3} =  {'8','16'};

subjects{4} = '12023_001';
niftiBaseName{4} = 'cm_12023_001_';
psirNiftiBaseName{4} = 'cm_12023_001';
T2star{4} = '14';
refScan(4) = 9; % scan before t2 structural
wholeheadPSIR{4} = '17';
distCorrectionRefSparse{4} = {'10','11'};
distCorrectionRefCont{4} = {'12','13'};
freeSurferName{4} = '12023_001';
sparseScans{4} =  {'8','15'};
contScans{4} =  {'9','16'};

subjects{5} = '11108_006';
niftiBaseName{5} = 'cm_11108_006_';
psirNiftiBaseName{5} = 'cm_11108_006';
T2star{5} = '17';
refScan(5) = 12; % scan before t2 structural
wholeheadPSIR{5} = '20';
distCorrectionRefSparse{5} = {'13','14'};
distCorrectionRefCont{5} = {'15','16'};
freeSurferName{5} = '11108_006';
sparseScans{5} =  {'8','18'};
contScans{5} =  {'12','19'};

subjects{6} = '11020_002';
niftiBaseName{6} = 'cm_11020_002_';
psirNiftiBaseName{6} = 'cm_11020_002';
T2star{6} = '14';
refScan(6) = 9; % scan before t2 structural
wholeheadPSIR{6} = '19';
distCorrectionRefSparse{6} = {'10','11'};
distCorrectionRefCont{6} = {'12','13'};
freeSurferName{6} = '11020_002';
sparseScans{6} =  {'08','15'};
contScans{6} =  {'09','18'};


cd([dataDir '/scanner/' subjects{iSubj}])

    %% check for filenames below 10 and add a zero before the number - make sure linking of stim files works later
    % str =
    % [token, remain] = strtok(str, ...)
scanFiles = dir;
for i = 1:length(scanFiles)
    str = scanFiles(i).name;
    strParts = strsplit(str,'_');
    checkScanNum = char(strParts(4));
    if 10>str2num(checkScanNum)
        
    end
end
    
  

system('ptoa -f -q -nii *.PAR')
% or try
!ptoa -f -q -nii *.PAR


if ~isempty(wholeheadPSIR{iSubj})
    % Move whole head PSIR
    mkdir(fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} wholeheadPSIR{iSubj} '*.nii']),fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    cd(fullfile(dataDir,'Anatomy/originals/',freeSurferName{iSubj}));
    system(sprintf('PSIR_script.sh . %s_%s_1 %s',psirNiftiBaseName{iSubj},wholeheadPSIR{iSubj},freeSurferName{iSubj}));
    %RUN RECON-AL_HIGHRES DIRECTLY IN TERMINAL
%     fprintf(['recon-all_highres ' freeSurferName{iSubj} ' ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' freeSurferName{iSubj} '_PSIR_pos_-.7_thr\n'])
    %     recon-all -all -s $SUBJECT -hires -i $IMAGE -expert $EXPERT_FILE
    system(sprintf('fslmaths %s_PSIR_pos_-.7_thr -mul 200 %s_PSIR_pos_-.7_thr_200',freeSurferName{iSubj},freeSurferName{iSubj}));
    % fprintf(['recon-all -all -s ' freeSurferName{iSubj} ' -hires -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' freeSurferName{iSubj} '_PSIR_pos_-.7_thr.nii' ' -expert ' fullfile(dataDir,'Anatomy','freesurfer','high_res_options.txt\n')])
    fprintf(['recon-all -all -s ' freeSurferName{iSubj} ' -hires -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' freeSurferName{iSubj} '_PSIR_pos_-.7_thr_200.nii' ' -expert ' fullfile(dataDir,'Anatomy','freesurfer','high_res_options.txt\n')])
    
end

correctDistortions = 1;
if correctDistortions
    mkdir(fullfile(dataDir,studyDir,subjects{iSubj}));
    cd(fullfile(dataDir,studyDir,subjects{iSubj}));
else
    mkdir(fullfile(dataDir,studyDir,subjects{iSubj},'_nonDistCorrected'));
    cd(fullfile(dataDir,studyDir,[subjects{iSubj} '_nonDistCorrected']));
end


mkdir('Etc')

mkdir('Anatomy')
mkdir('Raw')
mkdir('Raw/TSeries')
if correctDistortions
    mkdir('Distorted')
end
% mkdir ('T2_star')


%% Move in-plane T2 star
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} T2star{iSubj} '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Anatomy'))

%% Move the functional scans
if correctDistortions
    movefile(fullfile(dataDir,'scanner',subjects{iSubj},'*.nii'),fullfile(dataDir,studyDir,subjects{iSubj},'Distorted'))
    
    %% if distortion correction TEs are not interleaved
    cd('Distorted')

    
    check = 1;
    for n = 1:2
        if n == 1
            distCorrectionRef = distCorrectionRefSparse{iSubj};
        elseif n ==2
            distCorrectionRef = distCorrectionRefCont{iSubj};
        end
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
    end
    
    
    %% now do distortion correction
    % keyboard
    % opt = GetDistortionCorrectionParams(subjectID,distCorrectionScanNum,scanNum,refB0vol)
    optSparseFirst = GetDistortionCorrectionParams(subjects{iSubj},niftiBaseName{iSubj},'Sparse',distCorrectionRefSparse{iSubj},sparseScans{iSubj}{1, 1},'last');
    optSparseLast = GetDistortionCorrectionParams(subjects{iSubj},niftiBaseName{iSubj},'Sparse',distCorrectionRefSparse{iSubj},sparseScans{iSubj}{1, 2},'first');
    undistortProgram(optSparseFirst);
    undistortProgram(optSparseLast);
    optContFirst = GetDistortionCorrectionParams(subjects{iSubj},niftiBaseName{iSubj},'Cont',distCorrectionRefCont{iSubj},contScans{iSubj}{1, 1},'last');
    optContLast = GetDistortionCorrectionParams(subjects{iSubj},niftiBaseName{iSubj},'Cont',distCorrectionRefCont{iSubj},contScans{iSubj}{1, 2},'first');
    undistortProgram(optContFirst);
    undistortProgram(optContLast);
    cd('..')
    
    %% copy distortion corrected tseries to raw time series folder
    system('cp Distorted/dynB0map/*dynMod_U.nii Raw/TSeries/');
end

% crop last frame of reference EPI
!mkdir FNIRT
cd Raw/TSeries/


system(['fslroi ' niftiBaseName{iSubj} num2str(refScan(iSubj)) '_1_modulus_dynMod_U.nii lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);

!mv lastFrameEPI.nii ../../FNIRT
cd ../../

% crop T2*

cd Anatomy/
T2starFile = [niftiBaseName{iSubj} T2star{iSubj} '_1_modulus'];
system(sprintf('fslroi %s %s_crop 18 348 18 348 2 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));
% system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_stripped -f .2 -Z'])
system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_stripped -f .2 -Z'])

cd ../

%run FNIRT
cd FNIRT
system(sprintf('fnirtEpi2T2star lastFrameEPI %s_crop -separateFLIRT',T2starFile));

mrAlign
keyboard
% align cropped hi res t2* to PSIR .7_thr
% manually align the go to
% Compute Alignment->Advanced Alignment Menu
% then select Reverse Contrast (T2*)
% Click 'Compute Coarse Aligment' then
% Click 'Compute fine Aligment'

% once done
% save alignment to file and files (structural T2* images of the SAME size)

%% Below no longer needed - use function below
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
for i = 1:length(functionalNames)
    functionalNames{i} = ['../Raw/TSeries/' functionalNames{i}];
end
alignFunctional2HighResT2Star(sprintf('lastFrameEPI2%s_crop_resampled.affmat',T2starFile),[T2starFile '_crop.nii'],functionalNames);

cd ../


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
% mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[336 336 336]'); % this will depend on the resolution of the PSIR
mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[240 240 175]');

%apply FNIRT warping coefficient to surfaces
cd(fullfile(dataDir,studyDir,subjects{iSubj}));
fslApplyWarpSurfOFF(fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/',[niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_resampled_2_lastFrameEPI_warpcoef.nii']),...
    fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/','lastFrameEPI.nii'),...
    fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj},'surfRelax',[freeSurferName{iSubj} '_mprage_pp.nii']),...
    subjects{iSubj});

% add how to do freesurfer stuff
