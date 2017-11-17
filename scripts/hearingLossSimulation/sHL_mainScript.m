%% sHL_mainScript - scripted analysis for Simulated Hearing Loss (sHL) study 2017
% simulate hearing loss

%% TO DO

% check if current group is what you want to change to before changing

%% run this function when returning for analysis?
iSub = 1;
thisView = getMLRView;
[stimInfo, glmInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

% define and get subject info
% either cell for each subject or return single subject info - prob the later
subjectInfo = get_SubjectInfo_sHL(iSub);
% Subject ID
% flatmap names

% setup study direction (once)
% create folders needed
sHL_createStudyDirectory

% import, convert and move data
sHL_organiseData(studyDirectory,subject)
% per subject

preprocessData = sHL_preprocess;
% distortion correct
% linear alignment
% non-linear alignment

[thisView, concatedate] = script_setupmrLoadRet(thisView,groupNames);
% initiate mrLoadRet
% motion corerection
% group data

[thisView, glmData] = script_glmAnalysis(thisView);
% HRF = double gamma and box car
% All stims and 8 bins
% don't need to perform weighted mean on individual runs

thisView = script_importAnatomy(thisView);
% load in:
% reference EPI
% High resolution in-plane T2*
% surfaces
% create flatmaps
% MAKE ORIGINAL FLAT MAPS (LEFT AND RIGHT) USING MAKEFLAT AND NAME THEM
% [freeSurferName{iSubj} '_left_Flat.off'] AND [freeSurferName{iSubj} '_right_Flat.off']

[thisView, pRFdata] = script_pRFAnalysis(thisView);

%% CHANGE FLATMAP NAMES TO LEFT AND RIGHT

thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);
% create ROIs with the names:

% update this view
thisView = getMLRView;

% get condition names
conditionNames = cell(1,length(glmInfo.nStim));
for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
    analysisName = glmInfo.analysisNames_Scans{iAnal};
    conditionNames{iAnal} = get_analysisConditionNames(thisView,analysisName,'MotionComp',1);
end

%% Create subject data storage structure

%% Get analysis data from scans
% pass on to ROI analysis to restrict

%% get data from scans and groups
% save so don't need to load again
% change name to get_analysisData_GLM and save to analysisData.glm
scanData = getScanData_GLM(thisView,glmInfo.analysisNames_Scans,glmInfo.analysisNames_Groups,glmInfo.groupNames);

%% get data from ROIs
% save so don't need to load again
% Set group outside of script
roiData = script_getROIdata(thisView,scanData.scan_GLMdata,glmInfo.analysisBaseNames_Scans,Info.ROInames,glmInfo.analysisScanNum,'GLM');

%% Convert data to flatmap space

% export scan data
for iScan = 1:glmInfo.nScans
    for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
            analysisName = [glmInfo.analysisBaseNames_Scans{iAnal}, '_Scan_' mat2str(iScan)];
        for iSide = 1:length(subjectInfo.flatmapNames)
            thisView = script_covertData2FlatmapSpace(thisView,'MotionComp',analysisName,iScan,[],subjectInfo.flatmapNames{iSide});
        end
    end
end

%% export data for concat groups
% use getOverlay to get overlaynumber for analysis we want
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},glmInfo.analysisNames_Groups{iAnal},[],[41, 42],subjectInfo.flatmapNames{iSide});
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
        overlayNames= ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,[]);'])
        end
    end
end




%% restrict by roi
%% GROUPS - get data from ROIs
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




% roiAnalysis = script_ROIAnalysis(data_rightROI,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays');

%% save data
% saveLocation 
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
saveName = [subjectInfo.subjectID '_data.mat'];
% navigate to correct directory
% save roi analysis and data seperately due to size
save(saveName,'data','-v7.3','-nocompression')

left_ROIdata = data.Left.roiAnalysis;
right_ROIdata = data.Right.roiAnalysis;
saveName = [subjectInfo.subjectID '_ROIdata.mat'];
save(saveName,'left_ROIdata','right_ROIdata','-v7.3','-nocompression');
% save data to disk
% load later for group analysis


script_ROIAnalysis(left_ROIdata,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays',{'LeftAC'});

script_ROIAnalysis(right_ROIdata,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays',{'RightAC'});


plot_dbSLvsBetaWeight( left_ROIdata.LeftAC.roiAnalysis_glm_hrfDoubleGamma_nCons_32.roi_av_ratio , stimInfo.stimLevel_SL_mv , 50 )
plot_dbSLvsBetaWeight( right_ROIdata.RightAC.roiAnalysis_glm_hrfDoubleGamma_nCons_32.roi_av_ratio , stimInfo.stimLevel_SL_mv , 50 )
%% quit mrLoadRet
mrQuit()


script_GroupAnalysis(subjects)
% load subject ROI analysis from disk




