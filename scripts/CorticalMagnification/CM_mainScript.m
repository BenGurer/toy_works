%% CM_mainScript - scripted analysis for Cortical Magnification (CM) study 2017
% Measuring cortical magnification analysis
%
% load, sort and pre-process data
% GLM analysis
% pRF analysis
% Gradient reversals and ROI definition
% Cortical distance
% export data from volume to flat map space
% ROI analysis

% what are we interested in?
% HRF
% 

%% define subject
iSub = 1;

%% Get study parameters
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

%% Get subject info
subjectInfo = get_SubjectInfo_CM(iSub);
% Subject ID, flatmap names

%% if returning - move to folder and open mrLoadRet, and get/update thisView
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
mrLoadRet
thisView = getMLRView;

%% Pre-processing
CM_organiseData(studyDirectory,subject)
% import, convert and move data
% per subject
CM_preprocess;
% distortion correct
% linear alignment
% non-linear alignment
[thisView, concatedate] = script_setupmrLoadRet(thisView,groupNames);
% initiate mrLoadRet
% motion corerection
% group data

%% GLM analysis
[thisView, glmData] = script_glmAnalysis(thisView,glmInfo);
% HRF = double gamma and box car
% All stims and 8 bins
% don't need to perform weighted mean on individual runs

%% Load anatomy 
thisView = script_importAnatomy(thisView);
% load in:
% reference EPI
% High resolution in-plane T2*
% Surfaces
% create flatmaps
% MAKE ORIGINAL FLAT MAPS (LEFT AND RIGHT) USING MAKEFLAT AND NAME THEM
% [freeSurferName{iSubj} '_left_Flat.off'] AND [freeSurferName{iSubj} '_right_Flat.off']

%% pRF analysis
[thisView, pRFdata] = script_pRFAnalysis(thisView,pRFInfo);
% get info from glm to inform pRF 
% BOLD change between conditions
% average tuning curve sigma

%% Flatmap analysis
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);
% gradient reversals
% ROI defindition
% create ROIs with the names:

%% Get condition names
conditionNames = cell(1,length(glmInfo.nStim));
for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
    analysisName = glmInfo.analysisNames_Scans{iAnal};
    conditionNames{iAnal} = get_analysisConditionNames(thisView,analysisName,'MotionComp',1);
end

%% Convert data to flatmap space and average over depth
% export data from individual scans
for iScan = 1:glmInfo.nScans
    for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
            analysisName = [glmInfo.analysisBaseNames_Scans{iAnal}, '_Scan_' mat2str(iScan)];
        for iSide = 1:length(subjectInfo.flatmapNames)
            thisView = script_covertData2FlatmapSpace(thisView,'MotionComp',analysisName,iScan,[],subjectInfo.flatmapNames{iSide});
        end
    end
end

% export data from concatenated groups
% currently hard code overlay numbers - add to getStudyParams
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(subjectInfo.flatmapNames)        
        baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
        thisView = viewSet(thisView,'currentbase',baseNum);
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},glmInfo.analysisNames_Groups{iAnal},[],[1:33, 41, 42],subjectInfo.flatmapNames{iSide});
        end
    end
end

% average overdepth
for iSide = 1:length(subjectInfo.flatmapNames)
    thisView = script_averageAcrossDepths(thisView,[],[subjectInfo.flatmapNames{iSide}, 'Volume']);
end


%% GET DATA
%% SCANS - get data from Overlays
% create names to get data from overlays and save using structure side.Group.anal.data{iScan}
% get beta data from concat overlays to perform ROI anaylsis
q = char(39);
for iScan = 1:glmInfo.nScans
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
            eval(['data.' Info.Sides{iSide}, '.scans.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',conditionNames{iAnal},iScan);'])
        end
    end
end

%% GROUPS - get data from Overlays
conNamesString = [];
for iCon =1:length(conditionNames{1})
    if iCon == length(conditionNames{1})
    conNamesString = [conNamesString, conditionNames{1}{iCon}];
    else
    conNamesString = [conNamesString, conditionNames{1}{iCon}, ','];
    end
end

for iGroup = 1:length(glmInfo.groupNames)
    %         overlayNames= ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
    overlayNames = [];
    for iCon =1:length(conditionNames{1})
        overlayNames{iCon} = ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (' conditionNames{1}{iCon} ',0))'];
    end
    for iSide = 1:length(subjectInfo.flatmapNames)
        baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
        thisView = viewSet(thisView,'currentbase',baseNum);
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,[]);'])
        end
    end
end

%% restrict by roi
% GROUPS - get data from ROIs
% get data from left flatmap (set: left roi, left group, left base) save this data to left struct
% loop - save (side), selet:roi (side) and base (side), get data from: side.scans.anal.overlayData.data
% use eval
% save so don't need to load again
% Set group outside of script
% data_rightROI = script_getROIdata(thisView,Right,glmInfo.analysisBaseNames_Scans,{'FlatRightAC'},glmInfo.analysisScanNum,'overlays');
for iGroup = 1:length(glmInfo.groupNames) 
 for iSide = 1:length(Info.Sides)   
    eval(['dataVar = data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} ';']);
    eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
    eval(['ROI_data_' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} ' = script_getROIdata(thisView,dataVar,glmInfo.analysisNames_Groups,roiNames,[],' q 'overlays' q ');']);
 end
end

%% SCANS - get data from ROIs
for iSide = 1:length(Info.Sides)
    eval(['dataVar = data.' Info.Sides{iSide} ';']);
    eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
    eval(['ROI_data_' Info.Sides{iSide} ' = script_getROIdata(thisView,dataVar,glmInfo.analysisBaseNames_Scans,roiNames,glmInfo.analysisScanNum,' q 'overlays' q ');']);
end

%% perform ROI analysis
% save so don't need to load again
% compare binning for glm to averaging betas
% roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');
for iSide = 1:length(Info.Sides)
    eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
    eval(['data.' Info.Sides{iSide} '.roiAnalysis = script_ROIAnalysis(ROI_data_' Info.Sides{iSide} ',glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,' q 'overlays' q ',roiNames);']);
end

%% save data
% move to save location 
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));

% save data to structure with subject ID (unique name when loading for group analysis)
eval(['data_' subjectInfo.subjectID '.roi.left = data.Left.roiAnalysis;'])
eval(['data_' subjectInfo.subjectID '.roi.right = data.Right.roiAnalysis;'])

% save data
saveName = [subjectInfo.subjectID '_data.mat'];
% navigate to correct directory
% save roi analysis and data seperately due to size
eval(['save(saveName,data_' subjectInfo.subjectID ',' q '-v7.3' q ',' q '-nocompression' q ')']);

eval('roiData.left = data.Left.roiAnalysis;')
eval([subjectInfo.subjectID 'roiData.right = data.Right.roiAnalysis;'])
right_ROIdata = data.Right.roiAnalysis;
saveName = [subjectInfo.subjectID '_ROIdata.mat'];
save(saveName,'left_ROIdata','right_ROIdata','-v7.3','-nocompression');
% save data to disk
% load later for group analysis


cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
load([subjectInfo.subjectID '_ROIdata.mat']);
script_ROIAnalysis(left_ROIdata,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays',{'LeftAC'});

script_ROIAnalysis(right_ROIdata,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays',{'RightAC'});


plot_dbSLvsBetaWeight( left_ROIdata.LeftAC.roiAnalysis_glm_hrfDoubleGamma_nCons_32.roi_av_ratio , stimInfo.stimLevel_SL_mv , 50 )
plot_dbSLvsBetaWeight( right_ROIdata.RightAC.roiAnalysis_glm_hrfDoubleGamma_nCons_32.roi_av_ratio , stimInfo.stimLevel_SL_mv , 50 )
%% quit mrLoadRet
mrQuit()


script_GroupAnalysis(subjects)
% load subject ROI analysis from disk




