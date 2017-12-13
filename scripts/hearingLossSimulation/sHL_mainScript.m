%% sHL_mainScript - scripted analysis for Simulated Hearing Loss (sHL) study 2017
% simulate hearing loss

%% TO DO
% check if current group is what you want to change to before changing
% optimise prefit for pRF analysis
% - logical pCF and pTW range and suffient sampling of it
% - plus whats going on with non responsive voxels - ie pRF outside of stimulus range
% set up pRf info
% % add info to start of functions
% change gradient reversals to work with pRF analysis


iSub = 1;
q = char(39);

%% get study parameters
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

%% get subject info
% returns single subject info as a structure
subjectInfo = get_SubjectInfo_sHL(iSub);
% Subject ID
% flatmap names

%% move to subject folder, delete any current views, open mrLoadRet and get its view
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
deleteView(thisView);
mrLoadRet
thisView = getMLRView;
refreshMLRDisplay(thisView.viewNum);

%% organise subject data
sHL_organiseData(studyDirectory,subject);
% import, convert and move subject data

%% pre-process
sHL_preprocess;
% distortion correct
% linear alignment
% non-linear alignment

%% Setup mrLoadRet
[thisView, concatedate] = script_setupmrLoadRet(thisView,groupNames);
% initiate mrLoadRet
% motion corerection
% group data

%% Import anatomy
thisView = script_importAnatomy(thisView);
% load in:
% reference EPI
% High resolution in-plane T2*
% surfaces
% create flatmaps
% MAKE ORIGINAL FLAT MAPS (LEFT AND RIGHT) USING MAKEFLAT AND NAME THEM
% [freeSurferName{iSubj} '_left_Flat.off'] AND [freeSurferName{iSubj} '_right_Flat.off']

%% GLM analysis
thisView = script_glmAnalysis(thisView,glmInfo);
% HRF = double gamma and box car
% All stims and 8 bins
% don't need to perform weighted mean on individual runs

%% GLM grandient reversals
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);
% create ROIs with the names: GLMLeft, GLMright

%% Convert data to flatmap space and average over cortical depth
% export scan data
for iScan = 1:glmInfo.nScans
    for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
        analysisName = [glmInfo.analysisBaseNames_Scans{iAnal}, '_Scan_' mat2str(iScan)];
        for iSide = 1:length(subjectInfo.flatmapNames)
            thisView = script_covertData2FlatmapSpace(thisView,'MotionComp',analysisName,iScan,[],subjectInfo.flatmapNames{iSide});
        end
    end
end

% export group data
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(subjectInfo.flatmapNames)
        baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
        thisView = viewSet(thisView,'currentbase',baseNum);
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},glmInfo.analysisNames_Groups{iAnal},[],[1:33, 41, 42],subjectInfo.flatmapNames{iSide});
        end
    end
end

% average over cortical depth
for iSide = 1:length(subjectInfo.flatmapNames)
    thisView = script_averageAcrossDepths(thisView,[],[subjectInfo.flatmapNames{iSide}, 'Volume']);
end

%% get GLM data
% data has now be converted to flatmap space and averaged across cortical depth,
% the follow gets data from flatmap group, in flatmap space, from one layer that is the average across cortical depth.
% create names to get data from overlays and save using structure side.Group.anal.data{iScan}

% get condition names
conditionNames = cell(1,length(glmInfo.nStim));
for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
    analysisName = glmInfo.analysisNames_Scans{iAnal};
    conditionNames{iAnal} = get_analysisConditionNames(thisView,analysisName,'MotionComp',1);
end

% get data from SCANS
for iScan = 1:glmInfo.nScans
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
            eval(['data.' Info.Sides{iSide}, '.scanData.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',conditionNames{iAnal},iScan);'])
        end
    end
end

% get data from GROUPs
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

%% Restrict data by ROIs
% Now we need to restrict the data by the ROIs

% loop - save (side), selet:roi (side) and base (side), get data from: side.scans.anal.overlayData.data
% i.e. get data from left flatmap (set: left roi, left group, left base) save this data to left struct

% GROUPS
for iSide = 1:length(Info.Sides)
    baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
    thisView = viewSet(thisView,'currentbase',baseNum);
    for iGroup = 1:length(glmInfo.groupNames)
        eval(['groupDataVar = data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} ';']);
        eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
        eval(['ROI_data_' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} ' = script_getROIdata(thisView,groupDataVar,glmInfo.analysisNames_Groups,roiNames,[],' q 'overlays' q ');']);
%     eval(['ROI_data_' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} ' = script_getROIdata(thisView,dataVar,glmInfo.analysisNames_Groups,roiNames,[],' q 'overlays' q ');']);
%   eval(['data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} ' = script_getROIdata(thisView,dataVar,glmInfo.analysisNames_Groups,roiNames,[],' q 'overlays' q ');']);

    end
% end

    eval(['scanDataVar = data.' Info.Sides{iSide} ';']);    
    eval(['ROI_data_' Info.Sides{iSide} '.scanData = script_getROIdata(thisView,scanDataVar,glmInfo.analysisBaseNames_Scans,roiNames,glmInfo.analysisScanNum,' q 'overlays' q ');']);

% SCANS - get data from ROIs
% for iSide = 1:length(Info.Sides)
%     eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
%  eval(['data.' Info.Sides{iSide} '.scanData = script_getROIdata(thisView,dataVar,glmInfo.analysisNames_Groups,roiNames,[],' q 'overlays' q ');']);

end

%% GLM ROI analysis
% results from GLM roi analysis inform the pRF analysis so perform it first
%% perform ROI analysis
% save so don't need to load again
% compare binning for glm to averaging betas
% roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');
for iSide = 1:length(Info.Sides)
    eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
    eval(['data.' Info.Sides{iSide} '.roiNames.roiAnalysis = script_ROIAnalysis(ROI_data_' Info.Sides{iSide} ',glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,' q 'overlays' q ',roiNames);']);
 end


%% now create pRF restrict ROI in flat space and project through depths,
% can I script this?
% go to each flat base vol, define large ROI around HG, project through depths (ROIs>transform>expandROI([1 1 6])(replace)),
roi = viewGet(thisView,'roi','leftpRFrestrict');
roi = viewGet(thisView,'roi','pRFrestrict');
roi = viewGet(thisView,'roi','rightpRFrestrictexpand');
outputRoi = convertFromFlatVolumeToBase(roi);
thisView = viewSet(thisView,'newROI',outputRoi);
% then convert back to vol space, THEN:  ROIs>combine>(choose both ROIS),action:Union, newName:pRFrestrict
% name = pRFrestrict

%% pRF analysis
% get info from glm analysis: BOLD ratio, roi av TW
glmInfo.m = 0.0174;
glmInfo.b = -0.1176;
[thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,'rightpRFrestrictexpandVOL',1);
% get info from glm to inform pRF
% BOLD change between conditions
% average tuning curve sigma

%% pRF grandient reversals
% doesn't work?
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);
% create ROIs with the names: pRFleft, pRFright

%% get pRF data
% list overlay names and get using script_getOverlayData
for iGroup = 1:length(glmInfo.groupNames)
    %         overlayNames= ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
    
    overlayNames = 'averageDepthVol(pRF(PrefCentreFreq,0))';
    
    for iSide = 1:length(subjectInfo.flatmapNames)
        baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
        thisView = viewSet(thisView,'currentbase',baseNum);
        eval(['data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.', pRFinfo.analysisNames_Groups{iGroup}{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,[]);'])
        %         for iAnal = 1:length(glmInfo.analysisNames_Groups)
        %             eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,[]);'])
        %         end
    end
end

%% Restrict pRF data by ROI
% difference number of analysis in each group
% use the same format as glm - iAnal and loop through for each group
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(Info.Sides)
        eval(['dataVar = data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.' pRFinfo.analysisNames_Groups{iGroup}{iAnal}, '.overlayData;']);
        eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
        eval(['ROI_data_' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} ' = script_getROIdata(thisView,dataVar,pRFanalysis,roiNames,[],' q 'overlays' q ');']);
    end
end

%% pRF ROI analysis
% create function or add to: voxel comparisions;
% pCF scatter, pCF distribution, pCF correlation
% r = corr2(A,B)


%% Save data







% roiAnalysis = script_ROIAnalysis(data_rightROI,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'overlays');

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


%% scrap code here:$$
%% get data from scans and groups
% save so don't need to load again
% change name to get_analysisData_GLM and save to analysisData.glm
scanData = getScanData_GLM(thisView,glmInfo.analysisNames_Scans,glmInfo.analysisNames_Groups,glmInfo.groupNames);

%% get data from ROIs
% save so don't need to load again
% Set group outside of script
roiData = script_getROIdata(thisView,scanData.scan_GLMdata,glmInfo.analysisBaseNames_Scans,Info.ROInames,glmInfo.analysisScanNum,'GLM');

% conNamesString = [];
% for iCon =1:length(conditionNames{1})
%     if iCon == length(conditionNames{1})
%     conNamesString = [conNamesString, conditionNames{1}{iCon}];
%     else
%     conNamesString = [conNamesString, conditionNames{1}{iCon}, ','];
%     end
% end





