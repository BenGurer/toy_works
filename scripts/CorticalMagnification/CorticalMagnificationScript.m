%% export files with the same name and different scan number
% make cell array with subject info and then reference that

%% RUN FREESURFER FIRST

iSubj = 1;

epiDims = [128 128 24 361]; % dims of contin

% add check from startup

if ispc
    dataDir = 'N:/data';
elseif isunix
    dataDir = '/home/beng/data';
end
studyDir = 'CorticalMagnification';

sides = {'left','right'};
Sides = {'Left','Right'};

% Subject info
subjects{1} = '03644_012';
niftiBaseName{1} = 'pRFpilot2_';
T2star{1} = '16';
refScan{1} = '15'; % scan before t2 structural
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
refScan{2} = '8'; % scan before t2 structural
wholeheadPSIR{2} = '16';
distCorrectionRefSparse{2} = {'9','10'};
distCorrectionRefCont{2} = {'11','12'};
freeSurferName{2} = '12013_001';
sparseScans{2} =  {'7','14'};
contScans{2} =  {'8','15'};

subjects{3} = '12022_001';
psirSubject{3} = '12022_001';
niftiBaseName{3} = 'cm_12022_001_';
psirNiftiBaseName{3} = 'cm_12022_001';
T2star{3} = '13';
refScan{3} = '08'; % scan before t2 structural
wholeheadPSIR{3} = '17';
distCorrectionRefSparse{3} = {'09','10'};
distCorrectionRefCont{3} = {'11','12'};
freeSurferName{3} = '12022_001';
sparseScans{3} =  {'7','15'};
contScans{3} =  {'8','16'};

subjects{4} = '12023_001';
psirSubject{4} = '12023_001';
niftiBaseName{4} = 'cm_12023_001_';
psirNiftiBaseName{4} = 'cm_12023_001';
T2star{4} = '14';
refScan{4} = '9'; % scan before t2 structural
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
refScan{5} = '12'; % scan before t2 structural
wholeheadPSIR{5} = '20';
distCorrectionRefSparse{5} = {'13','14'};
distCorrectionRefCont{5} = {'15','16'};
freeSurferName{5} = '11108_006';
sparseScans{5} =  {'08','18'};
contScans{5} =  {'12','19'};

subjects{6} = '11020_002';
niftiBaseName{6} = 'cm_11020_002_';
psirNiftiBaseName{6} = 'cm_11020_002';
T2star{6} = '14';
refScan{6} = '09'; % scan before t2 structural
wholeheadPSIR{6} = '19';
distCorrectionRefSparse{6} = {'10','11'};
distCorrectionRefCont{6} = {'12','13'};
freeSurferName{6} = '11020_002';
sparseScans{6} =  {'08','15'};
contScans{6} =  {'09','18'};

% Import raw scans from subject 9 post motion correction
% Import to one group
subjects{7} = '08773_007';
niftiBaseName{7} = 'cm_08773_007_';
psirNiftiBaseName{7} = 'cm_08773_007';
T2star{7} = '16';
refScan{7} = '10'; % scan before t2 structural
wholeheadPSIR{7} = '17';
distCorrectionRefSparse{7} = {'14','15'};
distCorrectionRefCont{7} = {'21','22'};
freeSurferName{7} = '08773_007';
sparseScans{7} =  {'10','23'};
contScans{7} =  {'20',[]};

% subject 9 is subject 7's repeat session of 3 functional scans
subjects{9} = '08773_008';
niftiBaseName{9} = 'cm_08773_008_';
T2star{9} = '11';
refScan{9} = '10'; % scan before t2 structural
distCorrectionRefSparse{9} = {'12','13'};
distCorrectionRefCont{9} = {'14','15'};
sparseScans{9} =  {'10',[]};
contScans{9} =  {'09','18'};

subjects{8} = '09933_005';
niftiBaseName{8} = 'cm_09933_005_';
psirNiftiBaseName{8} = 'cm_09933_005';
T2star{8} = '15';
refScan{8} = '14'; % scan before t2 structural
wholeheadPSIR{8} = '14_scan1';
distCorrectionRefSparse{8} = {'07','08'};
distCorrectionRefCont{8} = {'09','10'};
freeSurferName{8} = '09933_005';
sparseScans{8} =  {'07_scan1','13'};
contScans{8} =  {'06','14'};

cd([dataDir '/scanner/' subjects{iSubj}])

%%
cd(fullfile(dataDir,studyDir,subjects{iSubj}));
    
    
    %% check for filenames below 10 and add a zero before the number - make sure linking of stim files works later
    % str =
    % [token, remain] = strtok(str, ...)
scanFiles = dir;
for id = 1:length(scanFiles)
    str = scanFiles(id).name;
    strParts = strsplit(str,'_');
    if length(strParts) > 1
    if strcmp([char(strParts(2)) '_' char(strParts(3))],subjects{iSubj})
        checkScanNum = char(strParts(4));
        numStrParts = length(strParts);
        if 10>str2num(checkScanNum) && numel(checkScanNum)==1
            newName = [char(strParts(1)) '_' char(strParts(2)) '_' char(strParts(3)) '_' ['0' checkScanNum] '_' char(strParts(5))];
            movefile(scanFiles(id).name,newName);
        end
    end
    end
end

% system('ptoa -f -q -nii *.PAR')
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

%% Pre-process functional
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
% if scan is continuous
system(['fslroi ' niftiBaseName{iSubj} refScan{iSubj} '_1_modulus_dynMod_U.nii lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);
% if scan is sparse
system(['fslroi ' niftiBaseName{iSubj} refScan{iSubj} '_1_modulus_dynMod_U.nii lastFrameEPI 72  1'])

!mv lastFrameEPI.nii ../../FNIRT
cd ../../


%% Alignment
cd Anatomy/
mrAlign
% align cropped hi res t2* to PSIR .7_thr
% manually align the go to
% Compute Alignment->Advanced Alignment Menu
% then select Reverse Contrast (T2*)
% Click 'Compute Coarse Aligment' then
% Click 'Compute fine Aligment'

% crop T2*
% T2* should be aligned to the wholehead PSIR at this point
T2starFile = [niftiBaseName{iSubj} T2star{iSubj} '_1_modulus'];
system(sprintf('fslroi %s %s_crop 18 348 18 348 2 24 0 1',T2starFile,T2starFile)); %% cropping t2 inplane structural image to match size of functional data
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));
system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_stripped -f .2 -Z'])

mrAlign
% skull stripping removes the s-form matrix so using set alignment to
% identity and save (source/t2*_stripped to destination/t2*)

cd ../

%run FNIRT
cd FNIRT
system(sprintf('fnirtEpi2T2star lastFrameEPI %s_crop -separateFLIRT',T2starFile));

% align functional data to high resolution t2* inplane structural image
% Below no longer needed - use function below
% Reload source as destination
% load lastFrameEPI as source
% Click Manual Alignment-> Set Alignment to Identity
% Save alignment to files (all functional data, including distorted)

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

%% Set up mrLoadRet
[sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1'); % looks in Raw/Tseries to find out how many scans there are
sessionParams.subject = subjects{iSubj};
sessionParams.description = 'CorticalMagnification';
sessionParams.operator = 'bg';

groupParams.description([1,3]) = {'Sparse, Run 1','Sparse, Run 2'};
groupParams.description([2,4]) = {'Continuous, Run 1','Continuous, Run 2'};
nScans = length(groupParams.name);

mrInit(sessionParams,groupParams,'makeReadme=0');

%% motion correction
thisView = newView;
refScanNum = viewGet(thisView,'scannum',sprintf('%s%s_1_modulus_dynMod_U.nii',niftiBaseName{iSubj},refScan{iSubj}));
[thisView, motionCompParams] = motionComp(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(1:nScans)]);
motionCompParams.baseFrame='last';
motionCompParams.baseScan = refScanNum;
[thisView, motionCompParams] = motionComp(thisView,motionCompParams);

%% concatenation of sparse data
thisView = viewSet(thisView,'curGroup','MotionComp');
params = getConcatParams_withNewGroupName(thisView,'ConcatenationSparse','defaultParams=1',['scanList=' mat2str([1 3])]);
[thisView, concatParamsSparse] = concatTSeries(thisView,params);
% [thisView, concatParams] = concatTSeries(thisView,[],'defaultParams=1',['scanList=' mat2str([1 3])]);

%% concatenation of continuous data
thisView = viewSet(thisView,'curGroup','MotionComp');
params = getConcatParams_withNewGroupName(thisView,'ConcatenationCont','defaultParams=1',['scanList=' mat2str([2 4])]);
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

% import PSIR and PD
thisView = viewSet(thisView,'newGroup','PSIR');
thisView = viewSet(thisView,'curGroup','PSIR');
thisView = importTSeries(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj},'TFE_PSIR_pos_-1_thr.nii')]);
thisView = newAnalysis(thisView,'dummy');
thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj},'TFE_PSIR_pos_-1_thr.nii')]);
thisView = viewSet(thisView,'overlaycolorrange',[.5 1]);
thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj},'TFE_PD_smooth7.nii')]);
%load PSIR as anatomy
thisView = loadAnat(thisView,'TFE_PSIR_pos_-1_thr.nii',fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj}));

%delayS = 2.5
%durationS = 2.5
% stimDurationMode = 'From File'
% stimDuration = []
% supersamplingMode = 'Set value'
% DesignSupersampling = 3
% acquistionDelay = 0.75
% numberContrasts = 0
% numberFtests = 0
% outputEstimatesAsOverlays = 1

%% add Ftest

thisView = getMLRView;
%% Statistical Analyses
% Loop between groups - NH and sHL
% GLM analysis
    concatenationGroup = {'ConcatenationSparse', 'ConcatenationCont'};
    
    concatenationGroup = {'Concatenation Sparse', 'Concatenation Cont'};
    hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
    for i = 1:length(concatenationGroup)        
        for ii = 1:length(hrfModel)
            analysisName = ['glm_' hrfModel{ii}];
            thisView = viewSet(thisView,'curGroup',concatenationGroup{i});
            [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
            glmParams.hrfModel = hrfModel{ii};
            [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
            glmParams.saveName = analysisName;
            glmParams.hrfParams.description = hrfModel{ii};
            switch hrfModel{ii}
                case 'hrfBoxcar'
            glmParams.hrfParams.delayS =  2.5;
            glmParams.hrfParams.durationS = 2.5;
                case 'hrfDoubleGamma'
                    
            end
            glmParams.scanParams{1}.stimDurationMode = 'From file';
            glmParams.scanParams{1}.supersamplingMode =  'Set value';
            glmParams.scanParams{1}.designSupersampling = 3;
            glmParams.scanParams{1}.acquisitionDelay = .75;
            glmParams.computeTtests = 1;
            glmParams.numberContrasts  = 1;
            glmParams.contrasts = ones(1,32);
            glmParams.numberFtests  = 1;
            glmParams.fTestNames{1, 1} = 'fTest - all conditions';
            glmParams.restrictions{1, 1} = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;...
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
            glmParams.alphaContrastOverlay = 'Uncorrected';
            glmParams.parametricTests = 1;
            glmParams.fweAdjustment = 0;
            glmParams.fdrAdjustment = 0;
            glmParams.outputStatistic = 0;
            glmParams.numberContrasts = 0;
            glmParams.outputEstimatesAsOverlays = 1;
            [thisView, glmParams] = glmAnalysis(thisView,glmParams);
            
            %Tonotopy analysis
            % Index max
            [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
            params.combineFunction='indexMax';
            params.nOutputOverlays=2;
            [thisView,params] = combineTransformOverlays(thisView,params);
            curOverlay=viewGet(thisView,'curOverlay');
            thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);
            
            % Weighted mean and corrected weighted mean
            params.combineFunction='weightedMeanStd';
            params.nOutputOverlays=4;
            [thisView,params] = combineTransformOverlays(thisView,params);
            curOverlay=viewGet(thisView,'curOverlay');
            thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
            thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
            thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);
            
            %% save analysis
            saveAnalysis(thisView,analysisName)            
        end         
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GLM Reverse Correlation - Continuous data concatenated %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This needs to be error checked

thisView = viewSet(thisView,'curGroup','ConcatenationCont');

thisView = viewSet(thisView,'curGroup','Concatenation Cont');
refreshMLRDisplay(thisView);

[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
glmParams.hrfModel = 'hrfDeconvolution';
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.scanParams{1, 1}.preprocess  = 'binStimFreq';
glmParams.hrfParams.hdrlenS = 15;
glmParams.numberContrasts = 8;
glmParams.componentsToTest = [0 1 1 1 1 0 0 0 0 0];
glmParams.numberEVs = 8;
glmParams.computeTtests = 1;
glmParams.numberFtests  = 1;
glmParams.fTestNames{1, 1} = 'fTest - all conditions';
glmParams.restrictions{1, 1} = [1,0,0,0,0,0,0,0;0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1];
glmParams.alphaContrastOverlay = 'Uncorrected';
glmParams.parametricTests = 1;
glmParams.fweAdjustment = 0;
glmParams.fdrAdjustment = 0;
glmParams.outputStatistic = 0;
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
glmParams.saveName = 'GLM_RevCorr_8bins';
glmParams.hrfParams.description = 'GLM Reverse Correlation - Cont Concat';
% glmParams.hrfParams.hdrlenS = 15;
% glmParams.numberContrasts = 8;
% glmParams.componentsToTest = [0 1 1 1 1 0 0 0 0 0];

[thisView, glmParams] = glmAnalysis(thisView,glmParams);

%% use interrogator getOverlayFromGlmAnalysis
% converts events to beta weights

%Tonotopy analysis
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);
% thisView = viewSet(thisView,'alphaOverlay',curOverlay,curOverlay-1);
% thisView = viewSet(thisView,'alphaOverlayExponent',0,curOverlay-1);
% thisView = viewSet(thisView,'overlaymin',1);

%% Use interrogator GLM plot to find average HRF shape of ROI

% import PSIR and PD
% Change to be aligned version
thisView = viewSet(thisView,'newGroup','PSIR');
thisView = viewSet(thisView,'curGroup','PSIR');
thisView = importTSeries(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj},[psirSubject{iSubj} '_PSIR_pos_-.7_thr.nii'])]);
thisView = newAnalysis(thisView,'dummy');
thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj},[psirSubject{iSubj} '_PSIR_pos_-.7_thr.nii'])]);
thisView = viewSet(thisView,'overlaycolorrange',[.5 1]);
thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj},[psirSubject{iSubj} '_PD_smooth7.nii'])]);
%load PSIR as anatomy
thisView = loadAnat(thisView,[psirSubject{iSubj} '_PSIR_pos_-.7_thr.nii'],fullfile(dataDir,'Anatomy/originals/',psirSubject{iSubj}));

%load reference EPI as anatomy
thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT'));
%load skull-stripped  EPI as overlay
[thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

% save('preProcessParams.mat','motionCompParams','concatParams');
save('preProcessParams.mat','motionCompParams');


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
thisView = getMLRView;

% THE REST REQUIRES THE FLAT MAPS
keyboard

% create flat maps if not done already
% repeat for distorted and non-distorted
params=[];
%import flat maps
params.path = fullfile(dataDir,'Anatomy/freesurfer/subjects/',freeSurferName{iSubj},'surfRelax');
for iSide=1:2
  params.anatFileName = [freeSurferName{iSubj} '_mprage_pp.nii'];
  params.flatRes=3;
  params.threshold = 1;
  params.baseValues = 'curvature';
  params.flatFileName = [freeSurferName{iSubj} '_' sides{iSide} '_WM_Flat_' flatName{iSubj,iSide} '.off'];
  params.outerCoordsFileName = [freeSurferName{iSubj} '_' sides{iSide} '_GM.off'];
  params.innerCoordsFileName = [freeSurferName{iSubj} '_' sides{iSide} '_WM.off'];
  params.curvFileName = [freeSurferName{iSubj} '_' sides{iSide} '_Curv.vff'];
  base = importFlatOFF(params);
  base.name = [freeSurferName{iSubj} '_' sides{iSide} '_Flat'];
  thisView = viewSet(thisView, 'newbase', base);
%   thisView = viewSet(thisView,'rotate',flatRotation(iSubj,iSide));
  thisView = viewSet(thisView,'corticalDepth',[0.2 0.8]);
end
refreshMLRDisplay(thisView);


flatmapInfo{5} = {'80_131_81_Rad60', 'dist'};
concatenationGroup = {'ConcatenationSparse', 'ConcatenationCont'};
functionalAnalysis = {'GLM_hrfDoubleGamma'};
mainOverlays(5) = 36;
ROInames = {'LeftPAC','RightPAC'};

%% Gradient reversals
% analysisType = viewGet(thisView,'analysisType');

% need to export to new group because over rights / need one for each
% flatmap/hemisphere
for iGroup = 1:length(concatenationGroup)
for iAnalysis = 1:length(functionalAnalysis)
for iSide=1:2
  % gradient reversals
  thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
  thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis{iAnalysis}));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['x' freeSurferName{iSubj} '_' sides{iSide} '_WM_Flat_' flatmapInfo{iSubj}{iSide}]));
  
refreshMLRDisplay(thisView);
  [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(mainOverlays(iSubj))]);
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

%% Split run analysis;
functionalAnalysis = {'GLM_BoxCar'};
ROInames = {'LeftPAC','RightPAC','PAC'};
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
glmParams.scanParams{iScan}.stimDurationMode = 'From File';
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
%% Get split run estiamtes
% need to loop over roiNum
analysisType = viewGet(thisView,'analysisType');
roiNum = 3;
for roiNum = 1:length(ROInames)
    roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
end
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
e = cell(1,4);
thisView = viewSet(thisView,'curgroup','MotionComp');
for iScan = 1:nScans
    thisView = viewSet(thisView,'curScan' ,iScan);
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',['GLM Box Car - Scan ' mat2str(iScan)]));
    analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',['GLM Box Car - Scan ' mat2str(iScan)]));
    glmData{iScan} = analysisData.d{iScan};
    analysisParams{iScan} = analysisData.params;
    r2data = analysisData.overlays(1).data{iScan};
    % get roi scan coords
    roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
    %get ROI estimates
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
    [e{iScan},volumeIndices] = getEstimates(glmData{iScan} ,analysisParams{iScan} ,volumeIndices');
    nVoxels = length(volumeIndices);
end
plotROIav_GLMBetaEstimates_SplitRuns_pRFfit(e{1},e{3});
plotROIav_GLMBetaEstimates_SplitRuns_pRFfit(e{2},e{4});
[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns(e{1},e{3});
[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns(e{2},e{4});


[ROIbetas, ROIStes] = plotROIav_GLMBetaEstimates_SplitRuns([e{1} e{3}],[e{2} e{4}]);

thisView = getMLRView;
nStim = 32;
analysisSaveName = {'GLM Box Car - ALL CONS - Scan ','pRF - ALL CONS - Scan '};
for iScan = 1:nScans    
thisView = viewSet(thisView,'curGroup','MotionComp',['curScan=' mat2str(iScan)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% GLM analysis %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
glmParams.hrfModel = 'hrfBoxcar';
glmParams.hrfParams.delayS =  2.5;
glmParams.hrfParams.durationS = 2.5;
glmParams.hrfParams.description = [analysisSaveName{1} mat2str(iScan)];
[thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
glmParams.scanParams{iScan}.stimDurationMode = 'From File';
glmParams.scanParams{iScan}.supersamplingMode =  'Set value';
glmParams.scanParams{iScan}.designSupersampling = 3;
glmParams.scanParams{iScan}.acquisitionDelay = .75;
glmParams.numberContrasts = 0;
glmParams.parametricTests = 0;
glmParams.outputEstimatesAsOverlays = 1; 
glmParams.saveName = [analysisSaveName{1} mat2str(iScan)];
[thisView, glmParams] = glmAnalysis(thisView,glmParams,['scanList=' mat2str(iScan)]);

%Tonotopy analysis
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:nStim+1])],['scanList=' mat2str(iScan)]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 nStim],curOverlay-1);

params.combineFunction='weightedMeanStd';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 nStim],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 nStim],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 nStim*1.25],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 nStim*1.25],curOverlay);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% pRF analysis %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
% pRFParams.saveName = [analysisSaveName{2} mat2str(iScan)];
% pRFParams.restrict = ['ROI: ' roi{1, roiNum}.name];
% pRFParams.pRFFit.supersampling = 1;
% pRFParams.pRFFit.fitHDR = 0;
% pRFParams.pRFFit.dispStimScan = iScan;
% [thisView, pRFParams] = pRF_auditory(thisView,pRFParams,['scanList=' mat2str(iScan)]);
end
% save analysis
for iScan = 1:nScans   
saveAnalysis(thisView,[analysisSaveName{1} mat2str(iScan)]);
% saveAnalysis(thisView,[analysisSaveName{2} mat2str(iScan)]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Get GLM Data %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thisView = getMLRView;
nScans = 4;
analysisSaveName = {'GLM Box Car - ALL CONS - Scan ','pRF - ALL CONS - Scan '};
ROInames = {'LeftPAC','RightPAC','PAC'};
roiNum = 3;
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});

for iScan = 1:nScans
    thisView = viewSet(thisView,'curGroup','MotionComp',['curScan=' mat2str(iScan)]);
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',[analysisSaveName{1} mat2str(iScan)]));
    analysisData{iScan} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName{1} mat2str(iScan)]));
    glmData{iScan} = analysisData{iScan}.d{iScan};
    analysisParams{iScan} = analysisData{iScan}.params;
    r2data = analysisData{iScan}.overlays(1).data{iScan};
    % get roi scan coords
    roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
    %get ROI estimates
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
    [estimate{iScan},volumeIndices] = getEstimates(glmData{iScan} ,analysisParams{iScan} ,volumeIndices');
    nVoxels = length(volumeIndices);
    % save data    
    splitData(iScan).glmData = glmData{iScan};
    splitData(iScan).estimates = estimate{iScan};
    splitData(iScan).betas = squeeze(estimate{iScan}.betas);
    splitData(iScan).betaSte = squeeze(estimate{iScan}.betaSte);
    splitData(iScan).r2 = permute(r2data(volumeIndices),[2 1]);
end
for iScan = 1:nScans
[val{iScan} index{iScan}] = max(splitData(iScan).betas);
end

concatenationGroup = {'ConcatenationSparse', 'ConcatenationCont'};

concatenationGroup = {'ConcatenationCont'};
hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
analysisSaveNameConcat = {'GLM_RevCorr_fTest','GLM_hrfDoubleGamma'};
for iGroup = 1:length(concatenationGroup)
    for i = 1:length(analysisSaveNameConcat)
        
        %% loop over concat group and get values
        %% get voxel restriction from concat of the same acq type
        %% use r2 or ftest for threshold
        %% take top 70% of voxels
        
        thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
        analysisData{i} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisSaveNameConcat{i}));
%         glmData{iGroup} = analysisData{iGroup}.d{i};
         glmData{iGroup} = analysisData{iGroup}.d{1};
        analysisParams{iGroup} = analysisData{iGroup}.params;
%         r2data = analysisData{iGroup}.overlays(1).data{i};
        r2data = analysisData{iGroup}.overlays(1).data{:};
        % get roi scan coords
        roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
        %get ROI estimates
        volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
        [e{iGroup},volumeIndices] = getEstimates(glmData{iGroup} ,analysisParams{iGroup} ,volumeIndices');
        nVoxels = length(volumeIndices);
        % save data
        concatData(iGroup).glmData = glmData{iGroup};
%         concatData(iGroup).glmData = 
        concatData(iGroup).estimates = e{iGroup};
        concatData(iGroup).betas = squeeze(e{iGroup}.betas);
        concatData(iGroup).betaSte = squeeze(e{iGroup}.betaSte);
        concatData(iGroup).r2 = permute(r2data(iGroup),[2 1]);
    end
end

fTestIndex = permute(analysisData{1, 1}.overlays(10).data{1}(volumeIndices),[2 1]);
fTestIndexThres = fTestIndex <= 0.05;
fTestIndexSum = sum(fTestIndexThres);
% fTestIndexThres = fTestIndex > (max(fTestIndex)-min(fTestIndex).*0.3);
% change to be percent of voxels?
r2IndexSparse = concatData(1).r2 > (max(concatData(1).r2)-min(concatData(1).r2)).*0.3;
r2MeanSparse = mean(concatData(1).r2);
nVoxelsSparse = sum(r2IndexSparse);
r2IndexContin = concatData(2).r2 > (max(concatData(2).r2)-min(concatData(2).r2)).*0.3;
r2MeanContin = mean(concatData(2).r2);
nVoxelsContin = sum(r2IndexContin);
figure;
scatter(index{2}(r2IndexContin),index{4}(r2IndexContin))

[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(splitData(1).estimates,splitData(3).estimates,r2IndexSparse);
[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(splitData(2).estimates,splitData(4).estimates,r2IndexContin);

[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(splitData(1).estimates,splitData(3).estimates,fTestIndexThres);
[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(splitData(2).estimates,splitData(4).estimates,fTestIndexThres);

plotROIav_GLMBetaEstimates_SplitRuns_pRFfit(splitData(1).estimates,splitData(3).estimates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Get pRF Data %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for iScan = 1:nScans
    analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName{1} mat2str(iScan)]));
    roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
    %get ROI estimates
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
    pCFdata = analysisData.overlays(3).data{iScan};
    pCFest{iScan} = pCFdata(volumeIndices);
    pTWdata = analysisData.overlays(4).data{iScan};
    pTWest{iScan} = pTWdata(volumeIndices);
end

thisView = getMLRView;

roiAC = viewGet(thisView,'roi','AC');
% roiAC = viewGet(thisView,'roi','RIGHT');
roiAC.scanCoords = getROICoordinates(thisView,roiAC);
analysisSaveName{1} = {'pRF_auditory_2','pRF_auditory_w_2'};
analysisSaveName{2} = {'pRF_auditory_fitHDR_2'};
% functionalAnalysis = {'GLM_BoxCar'};
%get ROI estimates
for iGroup = 1:length(concatenationGroup)
for i = 1:length(analysisSaveName{iGroup})  
thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
analysisData{i} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisSaveName{iGroup}{i}));
r2data = analysisData{i}.overlays(1).data{1};
volumeIndices = sub2ind(size(r2data),roiAC.scanCoords(1,:),roiAC.scanCoords(2,:),roiAC.scanCoords(3,:));
r2{i} = r2data(volumeIndices);
pCFdata = analysisData{i}.overlays(3).data{1};
pCFest{i} = pCFdata(volumeIndices);
pTWdata = analysisData{i}.overlays(4).data{1};
pTWest{i} = pTWdata(volumeIndices);
d{i} = analysisData{i}.d;

data(iGroup).d{i} = analysisData{i}.d;
data(iGroup).pCFest{i} = pCFdata(volumeIndices);
data(iGroup).pTWest{i} = pTWdata(volumeIndices);
data(iGroup).r2{i} = r2{i};
plotROIav_pRFEstimates(pCFest{i},pTWest{i},analysisData{i}.d{1}.scale(1,:),r2{i});
end
end
