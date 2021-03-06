%% Setup Subject info
iSubj = 2;

epiDims = [128 128 24 73]; % dims of functional scans


sides = {'left','right'};
Sides = {'Left','Right'};

% Subject info
% subjects{1} = '02344_034';
% niftiBaseName{1} = 'Sparse_02344_SPARSE_64dYN_1.5mmiso_TE40_SENSE_';
% wholeheadMPRAGE{1} = 'Sparse_02344_MPRAGE_SENSE_13_1';
% freeSurferName{1} = '02344_034';
% T2star{1} = 'Sparse_02344_High_res_t2__SENSE_12_1';
% refScan{1} = '11'; % scan before t2 structural
% flatmapName{1} = {'80_131_81_Rad60', '181_127_85_Rad60'};
subjects{1} = '02344_034';
niftiBaseName{1} = 'HL_02344_034_';
wholeheadMPRAGE{1} = '13';
freeSurferName{1} = '02344_034';
T2star{1} = '12';
refScan{1} = '11'; % scan before t2 structural
flatmapName{1} = {'80_131_81_Rad60', '181_127_85_Rad60'};
apScan{1} = 5;
paScan{1} = 6;
epiScans{1} = {'07', '09', '10', '11'};

subjects{2} = '12013_002';
niftiBaseName{2} = 'HL_12013_002_';
wholeheadMPRAGE{2} = '11';
freeSurferName{2} = '12013_002';
T2star{2} = '8';
refScan{2} = '05'; % scan before t2 structural
flatmapName{2} = {' ', ' '};
apScan{2} = 5;
paScan{2} = 6;
epiScans{2} = [03, 04, 08, 09];

subjects{3} = '12023_002';
niftiBaseName{3} = 'HL_12023_002_';
wholeheadMPRAGE{3} = '10';
freeSurferName{3} = '12023_002';
T2star{3} = '7';
refScan{3} = '04'; % scan before t2 structural
distCorrectionRef{3} = {'5','6'};
flatmapName{3} = {' ', ' '};

subjects{4} = '11108_007';
niftiBaseName{4} = 'HL_11108_007_';
wholeheadMPRAGE{4} = '10';
freeSurferName{4} = '11108_007';
T2star{4} = '7';
refScan{4} = '04'; % scan before t2 structural
flatmapName{4} = {' ', ' '};

subjects{5} = '13016_001';
niftiBaseName{5} = 'HL_113016_001_';
wholeheadMPRAGE{5} = '1';
freeSurferName{5} = '13016_001';
T2star{5} = '9';
refScan{5} = '02'; % scan before t2 structural
flatmapName{5} = {' ', ' '};

%% Move data from scanner to study/subject folders

cd([dataDir '/scanner/' subjects{iSubj}]) 

% convert from Par/Rec to niffti format
!ptoa -f -q -nii *.PAR

if ~isempty(wholeheadMPRAGE{iSubj})
    % Move whole head PSIR
    mkdir(fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} wholeheadMPRAGE{iSubj} '*.nii']),fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}));
    cd(fullfile(dataDir,'Anatomy/originals/',freeSurferName{iSubj}));        
    %RUN recon-all DIRECTLY IN TERMINAL 
    fprintf(['recon-all -subjid ' freeSurferName{iSubj} ' -i ' fullfile(dataDir,'Anatomy','originals',freeSurferName{iSubj}) '/' [niftiBaseName{iSubj} wholeheadMPRAGE{iSubj} '_1.nii'] ' -all']);

end

% make subject directory
mkdir(fullfile(dataDir,studyDir,subjects{iSubj}));
cd(fullfile(dataDir,studyDir,subjects{iSubj}));

mkdir('Etc')
mkdir('Anatomy')
mkdir('Raw')
mkdir('Raw/TSeries')
!mkdir FNIRT
!mkdir Distorted

% Move in-plane T2 star
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} T2star{iSubj} '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Anatomy'))
% movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Raw/TSeries'))
movefile(fullfile(dataDir,'scanner',subjects{iSubj},[niftiBaseName{iSubj} '*.nii']),fullfile(dataDir,studyDir,subjects{iSubj},'Distorted'))

% check for filenames below 10 and add a zero before the number - make sure linking of stim files works later
cd Raw/TSeries/
scanFiles = dir;
for id = 1:length(scanFiles)
    str = scanFiles(id).name;
    strParts = strsplit(str,'_');
    if length(strParts) > 1
    checkScanNum = char(strParts(4));
%         numStrParts = length(strParts);
        if 10>str2num(checkScanNum) && numel(checkScanNum)==1
            newName = [subjects{iSubj} '_0' char(strParts(4)) '.nii']
            movefile(scanFiles(id).name,newName);
        else
                        newName = [subjects{iSubj} '_' char(strParts(4)) '.nii']
            movefile(scanFiles(id).name,newName);
        end
   
    end
end


%% Pre-process
%% fsl topup Distortion Correction



% go to folder with files
% join fat shift epis

% copy files to folder distorted
% copy param file to distorted

system(['fslmerge -t topup_input.nii ' niftiBaseName{iSubj} num2str(apScan{iSubj}) '_1.nii ' niftiBaseName{iSubj} num2str(paScan{iSubj}) '_1.nii']);
system('topup --imain=topup_input.nii --datain=base_fat_shifted_AP_params.txt --subsamp=1 --fwhm=0  --out=topup_usable_output --iout=topup_check_output');

for iScan = 1:length(epiScans{iSubj})
    fileName = [niftiBaseName{iSubj} epiScans{iSubj}{iScan}];
    system(['applytopup --imain=' fileName '_1.nii --datain=base_fat_shifted_AP_params.txt --inindex=1 --topup=topup_usable_output --out=' fileName '_1_unwarped.nii --method=jac']);
end

system('cp Distorted/*unwarped.nii Raw/TSeries/');

system('cp Distorted/*_1.nii Raw/TSeries/');

% fslview output

% 
% system(['fslmerge -t topup_input.nii HL_12023_002_5_1.nii HL_12023_002_6_1.nii']);
           
% loop over scans

% 'applytopup --imain='niftiBaseName{iSubj} apScan{iSubj} '_1.nii --datain='<<same params file as above>> --inindex=1 --topup=<<same out name as above>> --out=<<output file for unwarped image>> --method=jac
% 
% 'applytopup --imain=<<image to be corrected>> --datain=<<same params file as above>> --inindex=1 --topup=<<same out name as above>> --out=<<output file for unwarped image>> --method=jac'
% 
% 'applytopup --imain=HL_12023_002_3_1.nii --datain=base_fat_shifted_params.txt --inindex=1 --topup=topup_usable_output --out=HL_12023_002_3_1_unwarped.nii --method=jac'


% crop last frame of reference EPI
system(['fslroi ' subjects{iSubj} '_' refScan{iSubj} '.nii lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);

system(['fslroi ' niftiBaseName{iSubj} refScan{iSubj} '_1_unwarped.nii lastFrameEPI ' num2str(epiDims(4)-1) ' 1']);

!mv lastFrameEPI.nii ../FNIRT
!mv lastFrameEPI.nii ../../FNIRT
cd ../../

% align T2* in-plane to T1 wholehead
% crop T2*
cd Anatomy/
mrAlign
keyboard
% if subject has been scanned before - first align MPRAGE to PSIR
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

% align EPI to T2* in-plane
% crop T2* in-plane to use with FNIRT and FLIRT
T2starFile = [niftiBaseName{iSubj} T2star{iSubj} '_1_modulus'];
T2starFile = [T2star{iSubj} '_modulus']; % subject 1
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
system(sprintf('cp %s_crop.nii ../FNIRT',T2starFile));
% system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_crop_stripped -f .2 -Z'])
system(['bet ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus ' niftiBaseName{iSubj} T2star{iSubj} '_1_modulus_stripped -f .2 -Z'])

mrAlign
keyboard
% skull stripping removes the s-form matrix so using set alignment to
% identity and save (t2*_stripped to t2*)

% Use FLIRT AND FNIRT to register functional to structural
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

%% Pre-processing DONE


%% Set up mrTools mrLoadRet
[sessionParams, groupParams] = mrInit([],[],'justGetParams=1','defaultParams=1'); % looks in Raw/Tseries to find out how many scans there are
sessionParams.subject = subjects{iSubj};
sessionParams.description = studyDir;
sessionParams.operator = 'bg';

groupParams.description([1,3]) = {'Hearing Loss Simulation, Run 1','Hearing Loss Simulation, Run 2'};
groupParams.description([2,4]) = {'Normal Hearing, Run 1','Normal Hearing, Run 2'};
nScans = length(groupParams.name);

mrInit(sessionParams,groupParams,'makeReadme=0');

% Motion correction

refScanNum = viewGet(thisView,'scannum',sprintf('%s%s.nii',niftiBaseName{iSubj},refScan{iSubj}));

thisView = newView;
refScanNum = viewGet(thisView,'scannum',sprintf('%s_%s.nii',subjects{iSubj},refScan{iSubj}));
[thisView, motionCompParams] = motionComp(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(1:nScans)]);
motionCompParams.baseFrame='last';
motionCompParams.baseScan = refScanNum;
[thisView, motionCompParams] = motionComp(thisView,motionCompParams);

% Concatenation of Hearing Loss Simulation data
thisView = viewSet(thisView,'curGroup','MotionComp');
params_ConcatenationHLsim = getConcatParams_withNewGroupName(thisView,'ConcatenationHLsim','defaultParams=1',['scanList=' mat2str([1 3])]);
[thisView, concatParamsSparse] = concatTSeries(thisView,params_ConcatenationHLsim);

% Concatenation of Normal Hearing data
thisView = viewSet(thisView,'curGroup','MotionComp');
params_ConcatenationNH = getConcatParams_withNewGroupName(thisView,'ConcatenationNH','defaultParams=1',['scanList=' mat2str([2 4])]);
[thisView, concatParamsCont] = concatTSeries(thisView,params_ConcatenationNH);

% link stim files to scans
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

% save('preProcessParams.mat','motionCompParams','concatParams');
save('preProcessParams.mat','motionCompParams','params_ConcatenationNH','params_ConcatenationHLsim');

%% Statistical Analyses
% Loop between groups - NH and sHL
% GLM analysis
    % Box car
    % Double gamma
% pRF analysis
    % Box car
    % Double gamma
% Loop over scans in motioncomp group
% GLM analysis
    % Box car
% pRF analysis
    % Box car

%% Get Data
% Get ROI
    % Define using anatomy
        % AC
    % Define using functional properties
        % gradient reversals
        % pTW    
% n = length(ROI)
% data.analysis.estimate(n)
% estimate = pCF,pTW,r2
% Threshold data
    % r2 values
    % f-test
    % beta weight
    
%% Visualise Data
% Within Condition
% Voxel pCF v pTW
% Binned pCF v averaged pTW

% Between Condition
% Condition A vs B
% Condition x: Run 1 vs 2
% Distrubtion of pCF estimates
% Binned pCF estimates
% Binned pCF v averaged pTW
% Voxel pCF agreement
% Voxel pTW agreement
% Voxel r2 agreement
% Voxel vector


thisView = getMLRView;
%% Statistical Analyses
concatenationGroup = {'ConcatenationHLsim', 'ConcatenationNH'};
ROIName = 'AC';
data = analysisScript_tonotopic(thisView,concatenationGroup,ROIName,nScans);

% Loop between groups - NH and sHL
% GLM analysis
    concatenationGroup = {'ConcatenationHLsim', 'ConcatenationNH'};
    hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
    for i = 1:length(concatenationGroup)        
        for ii = 1:length(hrfModel)
            analysisName = ['glm_' hrfModel{ii}];
            thisView = viewSet(thisView,'curGroup',concatenationGroup{i});
            [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
            glmParams.hrfModel = hrfModel{ii};
            [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
            glmParams.saveName = analysisName;
            glmParams.hrfParams.description = hrfModel;
            switch hrfmodel
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


% %load reference EPI as anatomy
thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT'));
%load skull-stripped  EPI as overlay
% [thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

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

fslApplyWarpSurfOFF(fullfile(dataDir,studyDir,subjects{iSubj},'FNIRT/',[T2star{iSubj} '_modulus_crop_resampled_2_lastFrameEPI_warpcoef.nii']),...
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

concatenationGroup = {'ConcatenationHLsim', 'ConcatenationNH'};
functionalAnalysis = ['glm_' hrfModel{iHRF}]
mainOverlays(1) = 39;
ROInames = {'LeftPAC','RightPAC'};

% analysisType = viewGet(thisView,'analysisType');
roiNum = 2;
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});


        analysisName = ['glm_' hrfModel{iHRF}]
% need to export to new group because over rights / need one for each
% flatmap/hemisphere
for iGroup = 1:length(concatenationGroup)
for iAnalysis = 1:length(functionalAnalysis)
    functionalAnalysis = ['glm_' hrfModel{iAnalysis}];
for iSide=1:2
      % gradient reversals
  thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
  thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['x' freeSurferName{iSubj} '_' sides{iSide} '_WM_Flat_' flatmapName{iSubj}{iSide}]));
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

%% make ROI's based on gradient reversals on flat maps - Label them LeftPAC and RightPAC

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
glmParams.fTestNames{1,iScan} = 'fTest 1';
glmParams.restrictions{1,iScan} = [1,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;...
    0,0,0,0,0,0,0,0;];
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
ROInames = {'LEFT','RIGHT','AC'};
roiNum = 3;
analysisType = viewGet(thisView,'analysisType');
% analysisParams = convertOldGlmParams(viewGet(thisView,'analysisParams'));
roiNum = 3;
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

'AC'
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
ROInames = {'LEFT','RIGHT','AC'};
roiNum = 3;
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
%%%%%%%%% Get Split data %%%%%%%
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

%%%%%%%%% Get Concat data %%%%%%%
hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
analysisSaveNameConcat = {'GLM_BoxCar'};
analysisData = cell(1,length(concatenationGroup));
glmData = cell(1,length(concatenationGroup));
analysisParams = cell(1,length(concatenationGroup));
r2data = cell(1,length(concatenationGroup));
e = cell(1,length(concatenationGroup));
for iGroup = 1:length(concatenationGroup)
    for i = 1:length(analysisSaveNameConcat)
        thisView = viewSet(thisView,'curgroup',concatenationGroup{iGroup});
        analysisData{iGroup} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisSaveNameConcat{i}));
        glmData{iGroup} = analysisData{iGroup}.d{1};
        analysisParams{iGroup} = analysisData{iGroup}.params;
        r2data = analysisData{iGroup}.overlays(1).data{1};
        % get roi scan coords
        roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
        %get ROI estimates
        volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
        [e{iGroup},volumeIndices] = getEstimates(glmData{iGroup} ,analysisParams{iGroup} ,volumeIndices');
        nVoxels = length(volumeIndices);
        % save data
        concatData(iGroup).glmData = glmData{iGroup};
        concatData(iGroup).estimates = e{iGroup};
        concatData(iGroup).betas = squeeze(e{iGroup}.betas);
        concatData(iGroup).betaSte = squeeze(e{iGroup}.betaSte);
        concatData(iGroup).r2 = permute(r2data(volumeIndices),[2 1]);
    end
end
r2Index = concatData(2).r2 > (max(concatData(2).r2)-min(concatData(2).r2)).*0.5;
r2MeanContin = mean(concatData(2).r2);
nVoxelsContin = sum(r2Index);

allIndex = true(1,length(concatData(2).betas ));

[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(splitData(1).estimates,splitData(3).estimates,r2Index);
[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage_pRFfit(splitData(2).estimates,splitData(4).estimates,r2Index);


thisView = getMLRView;
ROInames = {'LeftPAC','RIGHT'};

analysisType = viewGet(thisView,'analysisType');
% analysisParams = convertOldGlmParams(viewGet(thisView,'analysisParams'));
roiNum = 2;
% 
%     roiList = selectInList(thisView,'rois');
%     
%   roi = viewGet(thisView,'roi',iRoi);
%   
%   roi{roiNum} = viewGet(thisView,'roi','RightPAC');
  analysisType = viewGet(thisView,'analysisType');
roiNum = 2; 
for roiNum = 1:length(ROInames)
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
end
  roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
  
  analysisSaveName = 'pRF - ALL CONS - Scan ';
for iScan = 1:nScans
thisView = viewSet(thisView,'curGroup','MotionComp');
thisView = viewSet(thisView,'curScan', iScan);
[thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);

pRFParams.saveName = [analysisSaveName mat2str(iScan)];
pRFParams.restrict = ['ROI: ' roi{1, roiNum}.name];
pRFParams.pRFFit.supersampling = 1;
pRFParams.pRFFit.fitHDR = 0;
pRFParams.pRFFit.dispStimScan = iScan;
% [thisView, pRFParams] = pRF_auditory(thisView,[],'justGetParams=1',['scanList=' mat2str(iScan)]);
[thisView, pRFParams] = pRF_auditory(thisView,pRFParams,['scanList=' mat2str(iScan)]);
% save analysis
saveAnalysis(thisView,[analysisSaveName mat2str(iScan)])
% analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName mat2str(iScan)]));
end
for iScan = 1:nScans
analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName mat2str(iScan)]));

 roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
    %get ROI estimates
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
    pCFdata = analysisData.overlays(3).data{iScan};
    pCFest{iScan} = pCFdata(volumeIndices);
        pTWdata = analysisData.overlays(4).data{iScan};
    pTWest{iScan} = pTWdata(volumeIndices);
end
% 
%%% SCRAP CODE %%%%%%

HLpCFdata = pRF_auditory_HL.overlays(3).data{1, 1};
NHpCFdata = pRF_auditory_NH.overlays(3).data{1, 1};

HLpCF = HLpCFdata(volumeIndices);
NHpCF = NHpCFdata(volumeIndices);

figure; scatter(NHpCF,HLpCF)
e = cell(1,4);
roi{roiNum} = viewGet(thisView,'roi',ROInames{roiNum});
e = cell(1,4);
thisView = viewSet(thisView,'curgroup','MotionComp');
for iScan = 1:nScans
    thisView = viewSet(thisView,'curScan' ,iScan);
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',[analysisSaveName mat2str(iScan)]));
    analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName mat2str(iScan)]));
    glmData{iScan} = analysisData.d{iScan};
    analysisParams{iScan} = analysisData.params;
    r2data = analysisData.overlays(1).data{iScan};
    
    % get roi scan coords
    roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
    %get ROI estimates
    volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
    pCFdata = analysisData.overlays(36).data{iScan};
    pCFest{iScan} = pCFdata(volumeIndices);
    %     pCFdata = analysisData.overlays(38).data{iScan};
    %     pCFest{iScan} = pCFdata(volumeIndices);
    pTWdata = analysisData.overlays(37).data{iScan};
    pTWest{iScan} = pTWdata(volumeIndices);
    [e{iScan},volumeIndices] = getEstimates(glmData{iScan} ,analysisParams{iScan} ,volumeIndices');
    nVoxels = length(volumeIndices);
end
[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns(e{1},e{3});

[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns(e{2},e{4});
[ROIbetasSum{2}, ROIStesSum{2}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage(e{2},e{4});

% weighted by r2
figure;scatter(pCFest{2},pCFest{4});
xlim([0 32]); ylim([0 32]);
figure;scatter(pCFest{1},pCFest{3});

figure;scatter(pTWest{1},pTWest{3});

figure;scatter(pTWest{2},pTWest{4});

figure;scatter(pCFest{2},pTWest{2});
hold on
scatter(pCFest{4},pTWest{4});
xlim([0 32]); ylim([0 15]);
[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage(e{1},e{3});

[ROIbetasSum{1}, ROIStesSum{1}] = plotROIav_GLMBetaEstimates_SplitRuns_MovingAverage(e{2},e{4});

%% NH vs HL
thisView = getMLRView;

roiAC = viewGet(thisView,'roi','AC');
% roiAC = viewGet(thisView,'roi','RIGHT');
roiAC.scanCoords = getROICoordinates(thisView,roiAC);
concatenationGroup = {'ConcatenationHLsim', 'ConcatenationNH'};
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
% get voxel pCF,pTW,scaling
% threshold by r2? need a better measure of fit for non linear fitting
% or only include ones within stim range
% group by pCF range or by location? could be based on position on HG
% take mean and std of groups (pCF and pTW and Scaling)
% plot for NH and HL on same graph
% measure diff

% do split estimates
% plotROIav_pRFEstimates(pCFest{1},pTWest{1},analysisData{i}.d{1}.scale(1,:),r2{1});

% add difference between estimates as variable

r2Index = data(2).r2{1} > 0.2;
nVoxels = sum(r2Index);
data(iGroup).pCFest{i}(r2Index)

mean(data(iGroup).r2{i}(r2Index))
min(data(iGroup).r2{i}(r2Index))
max(data(iGroup).r2{i}(r2Index))
for i = 1:2
figure

VoxelStruct(i).conA = data(1).pCFest{i}(data(2).pCFest{1}>3 & data(2).pCFest{1}<33);
VoxelStruct(i).conB = data(2).pCFest{1}(data(2).pCFest{1}>3 & data(2).pCFest{1}<33);

g(i) = gramm('x',VoxelStruct(i).conA,'y',VoxelStruct(i).conB);
% g.geom_point('alpha',0.05)
g(i).geom_point()
g(i).stat_glm();

% g(i).stat_bin();
g(i).draw()
end

% nhIndex = data(2).pCFest{1}>3 & data(2).pCFest{1}<33;
nhIndex = length(data(2).pCFest{1});
for i = 1:2
figure
VoxelStruct(i).conB = data(1).pCFest{i}(nhIndex & r2Index);
VoxelStruct(i).conA = data(2).pCFest{1}(nhIndex & r2Index);
g(i) = gramm('x',VoxelStruct(i).conA,'y',VoxelStruct(i).conB);
% g.geom_point('alpha',0.05)
g(i).geom_point()
g(i).stat_glm();
 g(i).set_names('x','Condition A - Frequency (ERB)','y','Condition B - Frequency (ERB)')
% g(i).stat_bin();
g(i).draw()
end

unpacked.pCF = [data(1).pCFest{1}(r2Index) data(1).pCFest{2}(r2Index) data(2).pCFest{1}(r2Index)];
% unpacked.con = [repmat(1,1,length(data(1).pCFest{1})) repmat(2,1,length(data(1).pCFest{1})) repmat(3,1,length(data(1).pCFest{1}))];
unpacked.con = [repmat({'Condition B - pRF'},1,sum(r2Index)) repmat({'Condition B - pRF_mod'},1,sum(r2Index)) repmat({'Condition A'},1,sum(r2Index))];


figure
hist  = gramm('x',unpacked.pCF,'color',unpacked.con);
hist.stat_bin('nbins',8);
hist.set_names('x','Frequency (ERB)','y','Number of Voxels','color','Condition:')
% hist.geom_label()
hist.draw()

hist = [];
figure
hist  = gramm('x',unpacked.pCF,'color',unpacked.con);
hist.set_names('x','Frequency (ERB)','y','Frequency','color','Condition:')
hist.stat_density();
hist.draw()

% violin pCF vs pTW?

%%%%
% make stimulus plot into a gramm object
figure
x = 0:40;
mu = pCFest{1}(2000);
sigma = pTWest{1}(2000);
rfModel = exp(-(((x-mu).^2)/(2*(sigma^2))));