function opt = GetDistortionCorrectionParams(subjectID,BaseName,acqType,distCorrectionScanNums,scanNum,refB0vol)


% subjectID = '111'; % subjectID needs to be a string
% distCorrectionScanNums % Scan numbers of distortion correction scans Needs to be a cell of two fields contraining strings
% acqType % Acquistion type. Needs to be a string
% refB0vol % first or last. Needs to be a string
if ispc
    dataDir = 'N:/data';
elseif isunix
    dataDir = '/home/beng/data';
end
studyDir = '/CorticalMagnification/';


%% load template
load(fullfile(dataDir,studyDir,'DistCorrectionConfigTemplate.mat'))

%% Change require files

% Change directory
opt.dataPath = [dataDir studyDir subjectID '/Distorted/']; % save as char

% Choose mask 
opt.maskPath  = [dataDir studyDir subjectID '/Distorted/' BaseName distCorrectionScanNums{1} '_modulus_firstFrame_brain.nii'];

% select scan to distortion correct
opt.dyn.modPath = {[dataDir studyDir subjectID '/Distorted/' BaseName scanNum '_1_modulus.nii']}; % save as cell
opt.dyn.phasePath = {[dataDir studyDir subjectID '/Distorted/' BaseName scanNum  '_1_phase.nii']}; % save as cell


% Change distortion correct scan
opt.B0.modPath = [dataDir studyDir subjectID '/Distorted/' BaseName distCorrectionScanNums{1} '_' distCorrectionScanNums{2} '_modulus.nii']; % save as char
opt.B0.phasePath =[dataDir studyDir subjectID '/Distorted/' BaseName distCorrectionScanNums{1} '_' distCorrectionScanNums{2} '_phase.nii']; % save as char

% Change ref BO volume (closest scan in functional data to distortion correction scans)
opt.dyn.refB0vol = refB0vol; % save as string
opt.outputDir = [dataDir studyDir subjectID '/Distorted/' 'dynB0map/'];

opt.outname = [dataDir studyDir subjectID '/Distorted/' ['config' acqType refB0vol '.mat']];

save([dataDir studyDir subjectID '/Distorted/' ['config' acqType refB0vol '.mat']],'opt')