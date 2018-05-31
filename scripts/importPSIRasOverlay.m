function thisView = importPSIRasOverlay(thisView,Info,subjectInfo)
    %
    %   usage: importPSIRasOverlay(thisView,Info,subjectInfo)
    %      by: Ben Gurer
    %    date: 05/22/2018
    % purpose: cimport PSIR as an overlay and base anatomy
    %   input: mrView, Study information and subject information
    %  output: mrView with PSIR overlay and anotmay
    %
    
    % move to anatomy folder
    cd(fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName))
    % rename PSIR file so it works with matlab
    eval (['!cp ', [subjectInfo.freeSurferName, '_PSIR_pos_-.7_thr.nii'], ' ' [subjectInfo.freeSurferName, '_PSIR_mrTools.nii']])
    
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
    
    % import PSIR and PD
    thisView = viewSet(thisView,'newGroup','PSIR');
    thisView = viewSet(thisView,'curGroup','PSIR');
    thisView = importTSeries(thisView,[],'defaultParams=1',['pathname=' fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName,[subjectInfo.freeSurferName, '_PSIR_mrTools.nii'])]);
    thisView = newAnalysis(thisView,'dummy');
    thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName,[subjectInfo.freeSurferName, '_PD_smooth7.nii'])]);
    thisView = importOverlay(thisView,[],'defaultParams=1',['pathname=' fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName,[subjectInfo.freeSurferName, '_PSIR_mrTools.nii'])]);
    thisView = viewSet(thisView,'overlaycolorrange',[.5 1]);
    %load PSIR as anatomy
    thisView = loadAnat(thisView,[subjectInfo.freeSurferName, '_PSIR_mrTools.nii'],fullfile(Info.dataDir,'Anatomy/originals/',subjectInfo.freeSurferName));

end