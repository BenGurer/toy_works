function thisView = script_importAnatomy(thisView,Info,subjectInfo)


% %load reference EPI as anatomy
thisView = loadAnat(thisView,'lastFrameEPI.nii',fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'FNIRT'));
%load skull-stripped  EPI as overlay
% [thisView,params] = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(dataDir,studyDir,subjects{iSubj},'/FNIRT/lastFrameEPI_stripped.nii')]);

% script running glm on each scan individually and performing split
% analysis

%IMPORT  FREESURFER SURFACE
cd(fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName));
% mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[336 336 336]'); % this will depend on the resolution of the PSIR
mlrImportFreeSurfer('defaultParams=1','volumeCropSize=[240 240 175]');

%apply FNIRT warping coefficient to surfaces
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
fslApplyWarpSurfOFF(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'FNIRT/',[subjectInfo.niftiBaseName subjectInfo.T2star '_1_modulus_crop_resampled_2_lastFrameEPI_warpcoef.nii']),...
    fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID,'FNIRT/','lastFrameEPI.nii'),...
    fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName, 'surfRelax',[subjectInfo.freeSurferName '_mprage_pp.nii']),...
    subjectInfo.subjectID);

%import surfaces (this step has not been made scriptable yet) 

for iSide=1:2
  base = importSurfaceOFF(fullfile(Info.dataDir,'Anatomy/freesurfer/subjects/',subjectInfo.freeSurferName,'surfRelax',...
   [subjectInfo.freeSurferName '_' Info.sides{iSide} '_GM.off']));
  thisView = viewSet(thisView, 'newbase', base);
  thisView = viewSet(thisView,'corticalDepth',[0.2 0.8]);
end

% save view and quit
mrSaveView(thisView);
end