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
iSub = 5;
q = char(39);


%% Get study parameters
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

% use is pc to set data directory - could do in cm_setupStduyparams
% Info.dataDir
Info.dataDir = '/Volumes/DataDisk/data';

%% Get subject info
subjectInfo = get_SubjectInfo_CM(iSub);

% Subject ID, flatmap names
saveName = [subjectInfo.subjectID '_data.mat'];

%% move to subject folder, delete any current views, open mrLoadRet and get its view
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
% deleteView(thisView);

%% either load data or mrView
if exist(saveName, 'file')
    load(saveName)
else
    data = struct;
end

mrLoadRet
thisView = getMLRView;

refreshMLRDisplay(thisView.viewNum);

%% organise subject data
% sHL_organiseData(Info, subjectInfo);
% import, convert and move subject data

%% pre-process
% sHL_preprocess(Info, subjectInfo, 0);
% distortion correct
% linear alignment
% non-linear alignment

%% Setup mrLoadRet
% [thisView, concatedate] = script_setupmrLoadRet(Info,subjectInfo,glmInfo);
% initiate mrLoadRet
% motion corerection
% smooth (optional)
% group data
% concatenate runs

%% open View
mrLoadRet

%% Import anatomy
thisView = script_importAnatomy(thisView,Info,subjectInfo);
% load in:
% reference EPI
% High resolution in-plane T2*
% surfaces

% rotate surfaces and then get view
thisView = getMLRView;

% Import PSIR
thisView = importPSIRasOverlay(thisView,Info,subjectInfo);

% import flatmaps without FNIRT coords (after flatmaps have been created)
thisView = script_importFlatmaps(thisView,Info,subjectInfo);

%% Tonotopic analysis

%% GLM Analysis with Boxcar HRF
% Perform GLM analysis with boxcar model of HRF
% hrf = Box Car
% concatenated runs only
% [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)
thisView = script_glmAnalysis(thisView,glmInfo,glmInfo.groupNames,{'hrfBoxcar'},1,1);


%% Make flatmaps
% radius = 55
% centre on HG (use R2 and f-test to guide)
% use default names
% resolution = 3; method = mrFlatMesh
% rotate flatmaps for easy viewing (do before exporting to flatmap space)

%% ROI CREATION
% create ROIs with the names:
% LeftAR and RightAR: Group=Sparse; Analysis=glm_hrfboxcar; Overlay=f-test - set to 0.005
% Create ROI - continuous voxels;
% ROIs>transform>expandROI([3 3 3])(convolves ROI with a sphere with a diameter of 2^3 voxels)
% name=(Left or Right)ARexp
% Project through depths 0.3 to 0.7 to remove voxels outside of grey matter for ALL ROIs
% combine LeftAR and RightAR = AR & combine LeftARexp and RightARexp = ARexp

% get view so we have the ROIs
thisView = getMLRView;

%% HRF estimate
% estimate hrf using deconvolution
thisView = script_hrfAnalysis(thisView,glmInfo.groupNames{2});

%% get av HRF estimate for Auditory Responsive (AR) ROI
% use results for GLM
% save result to data
[ thisView, data.hrf.x_doubleGamma, data.hrf.x_Gamma, data.hrf.x_dGamma, data.hrf.deconv] = script_hrfROIAnalysis(thisView,'AR',glmInfo);
% save to analysis structures to pass to functions
glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;

% Estimate and fit average ROI tuning width using recentre method from deconvolution data
% [ thisView, x_g, x_r, tw_Deconv, estimate, ~, nVoxels] = script_centredTWROIAnalysis(thisView,'AR',glmInfo);

%% GLM Analysis with Double Gamma HRF
% Perform GLM analysis with double gamma model of HRF
% hrf = Box Car
% All runs
% [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)
thisView = script_glmAnalysis(thisView,glmInfo,glmInfo.groupNames,{'hrfDoubleGamma'},1,1);


%% GLM grandient reversals
% calculate gradient reversals using GLM (double gamma) analysis.
% pCF estimation = Juliensdebiased method
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase, glmInfo.analysisNames_Groups{2},'[18 18 21]');

%% ROI CREATION
% get f-test maps in to flatspace and average overdepth
for iSide = 1:length(subjectInfo.flatmapNames)
    groupName = 'ConcatenationSparse';
    if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
        thisView = viewSet(thisView,'curgroup',groupName);
    end
    analysisName = 'glm_hrfBoxcar';
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
    overlayNum = viewGet(thisView,'overlayNum','P [F (fTest - all conditions)]');
    thisView = script_covertData2FlatmapSpace(thisView,'ConcatenationSparse',analysisName,[],overlayNum,subjectInfo.flatmapNames{iSide});
end
% average over cortical depth
for iSide = 1:length(subjectInfo.flatmapNames)
    thisView = script_averageAcrossDepths(thisView,[],[subjectInfo.flatmapNames{iSide}, 'Volume']);
end

%% ROI CREATION
% create ROIs with the names:
% LeftGR, RightGR based on gradient reversals, unsmoothed tonotopic maps and f-test maps

%% Convert GLM data to flatmap space and average over cortical depth
% [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
% export scan data
for iScan = 1:glmInfo.nScans
    for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
        analysisName = [glmInfo.analysisBaseNames_Scans{iAnal}, '_Scan_' mat2str(iScan)];
        for iSide = 1:length(subjectInfo.flatmapNames)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.scanGroupName,analysisName,iScan,[],subjectInfo.flatmapNames{iSide});
        end
    end
end

% export group data
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(subjectInfo.flatmapNames)
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},glmInfo.analysisNames_Groups{iAnal},[],[],subjectInfo.flatmapNames{iSide});
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
getConditionNames = cell(1,length(glmInfo.nStim));
% save condition names
if isfield(data, 'conditions')
    conditionNames = data.conditions;
else
    
    analysisName = glmInfo.analysisNames_Scans{1};
    conditionNames{1} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
    analysisName = glmInfo.analysisNames_Scans{3};
    conditionNames{2} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
    data.conditions = conditionNames;
end

% % get data from SCANS
% analysisName = 'combineTransformOverlays';
% for iSide = 1:length(subjectInfo.flatmapNames)
%     groupName = [subjectInfo.flatmapNames{iSide} 'Volume'];
%     thisView = viewSet(thisView,'curgroup',groupName);
%     thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
%     for iScan = 1:glmInfo.nScans
%         for iAnal = 1:length(glmInfo.hrfModel)*length(glmInfo.nStim)
%             % define overlay names
%             % First, deteremine how many stimuli
%             if glmInfo.analysisNStim{iAnal} == length(conditionNames{1});
%                 overlayNames = cell(1,length(conditionNames{1}));
%                 for iCon = 1:length(conditionNames{1})
%                     overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisBaseNames_Scans{iAnal} '_Scan_' mat2str(iScan) ' (' conditionNames{1}{iCon} ',0))'];
%                 end
%             else
%                 overlayNames = cell(1,length(conditionNames{2}));
%                 for iCon = 1:length(conditionNames{2})
%                     overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisBaseNames_Scans{iAnal} '_Scan_' mat2str(iScan)  ' (' conditionNames{2}{iCon} ',0))'];
%                 end
%             end
%             overlayData = get_overlayData(thisView,overlayNames);
%             eval(['data.' Info.Sides{iSide}, '.scanData.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} =  overlayData']);
%         end
%     end
% end

analysisName = 'combineTransformOverlays';
for iSide = 1:length(subjectInfo.flatmapNames)
    groupName = [subjectInfo.flatmapNames{iSide} 'Volume'];
    thisView = viewSet(thisView,'curgroup',groupName);
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
    for iScan = 1:glmInfo.nScans
        for iAnal = 1:length(glmInfo.analysisNames_nCons)
            % define overlay names
            % First, deteremine how many stimuli
            if glmInfo.analysisNStim{iAnal} == length(conditionNames{1});
                overlayNames = cell(1,length(conditionNames{1}));
                for iCon = 1:length(conditionNames{1})
                    overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' conditionNames{1}{iCon} ',0))'];
                end
            else
                overlayNames = cell(1,length(conditionNames{2}));
                for iCon = 1:length(conditionNames{2})
                    overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan)  ' (' conditionNames{2}{iCon} ',0))'];
                end
            end
            overlayData = get_overlayData(thisView,overlayNames);
            eval(['data.' Info.Sides{iSide}, '.scans.', glmInfo.analysisNames_nCons{iAnal}, '.betas{iScan} =  overlayData']);
        end
    end
end

%% get data from GROUPs
% define vocel property overlay names - index max, centriod, spread, julienCentriod, julienTuningWidth
voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
for iSide = 1:length(subjectInfo.flatmapNames)
    baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
    thisView = viewSet(thisView,'currentbase',baseNum);
    for iGroup = 1:length(glmInfo.groupNames)
        for iAnal = 1:length(glmInfo.analysisNames)
            % define overlay names
            % betas overlay names
            overlayNames = cell(1,length(conditionNames{1}));
            conNamesString = [];
            for iCon =1:length(conditionNames{1})
                overlayNames{iCon} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' conditionNames{1}{iCon} ',0))'];
                if iCon == 1
                    conNamesString  = [conNamesString conditionNames{1}{iCon}];
                else
                    conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                end
            end
            
            % r2 overlay name
            r2OverlayName = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (r2,0))'];
            
            % voxel property estmate overlay names
            voxelPropertyOverlayName = cell(1,length(voxelPropertyNames));
            for iName = 1:length(voxelPropertyNames) - 1
                voxelPropertyOverlayName{iName} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (Ouput ' num2str(iName) ' - weightedMeanStd(' conNamesString '),0))'];
            end
            voxelPropertyOverlayName{end} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (Ouput 1 - indexMax(' conNamesString '),0))'];
            
            % get overlay data using get_overlayData - outputs a structure with fields: .data & .name
            % R2
            clear tempData
            tempData = get_overlayData(thisView,r2OverlayName);
            eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.r2 = tempData;']);
            
            % beta weights
            clear tempData
            tempData = get_overlayData(thisView,overlayNames);
            eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.betas = tempData;']);
            
            % voxel property estiamtes
            for iName = 1:length(voxelPropertyNames)
                clear tempData
                tempData = get_overlayData(thisView,voxelPropertyOverlayName{iName});
                eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.', voxelPropertyNames{iName},' = tempData;']);
            end
            
        end
    end
end

% %% get data from GROUPs
% % define vocel property overlay names - index max, centriod, spread, julienCentriod, julienTuningWidth
% voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
% for iSide = 1:length(subjectInfo.flatmapNames)
%     baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
%     thisView = viewSet(thisView,'currentbase',baseNum);
%     for iGroup = 1:length(glmInfo.groupNames)
%         for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
%             % define overlay names
%             % betas overlay names
%             overlayNames = cell(1,length(conditionNames{1}));
%             conNamesString = [];
%             for iCon =1:length(conditionNames{1})
%                 overlayNames{iCon} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames_Groups{iAnal} ' (' conditionNames{1}{iCon} ',0))'];
%                 if iCon == 1
%                     conNamesString  = [conNamesString conditionNames{1}{iCon}];
%                 else
%                     conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
%                 end
%             end
%             
%             % r2 overlay name
%             r2OverlayName = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames_Groups{iAnal} ' (r2,0))'];
%             
%             % voxel property estmate overlay names
%             voxelPropertyOverlayName = cell(1,length(voxelPropertyNames));
%             for iName = 1:length(voxelPropertyNames) - 1
%                 voxelPropertyOverlayName{iName} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames_Groups{iAnal} ' (Ouput ' num2str(iName) ' - weightedMeanStd(' conNamesString '),0))'];
%             end
%             voxelPropertyOverlayName{end} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames_Groups{iAnal} ' (Ouput 1 - indexMax(' conNamesString '),0))'];
%             
%             % get overlay data using get_overlayData - outputs a structure with fields: .data & .name
%             % R2
%             clear tempData
%             tempData = get_overlayData(thisView,r2OverlayName);
%             eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.r2 = tempData;']);
%             
%             % beta weights
%             clear tempData
%             tempData = get_overlayData(thisView,overlayNames);
%             eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.betas = tempData;']);
%             
%             % voxel property estiamtes
%             for iName = 1:length(voxelPropertyNames)
%                 clear tempData
%                 tempData = get_overlayData(thisView,voxelPropertyOverlayName{iName});
%                 eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.', voxelPropertyNames{iName},' = tempData;']);
%             end
%             
%         end
%     end
% end

%% Restrict data by ROIs
% Now we need to restrict the data by the ROIs
% first get view so we have the ROIS
thisView = getMLRView;

% loop/save (side), select:roi and base, get data from: side.scans.anal.overlayData.data
% i.e. get data from left flatmap (set: left roi, left group, left base) save this data to left struct

% GROUPS

q = char(39);
for iSide = 1:length(Info.Sides)
    eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
    baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
    thisView = viewSet(thisView,'currentbase',baseNum);
    
    for iROI = 1:length(ROInames)
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = struct;']);
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
        
        % Restrict Group data
        for iGroup = 1:length(glmInfo.groupNames)
            for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
                %                 eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.overlayData;']);
                eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
                analysisName = glmInfo.analysisNames_Groups{iAnal};
                eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas = script_getROIdata(thisView,groupDataVar,analysisName,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi,[],' q 'overlays' q ');']);
            end
        end
        
        % Restrict Scan data
        eval(['scanDataVar = data.' Info.Sides{iSide} '.scanData;']);
        %     eval(['ROI_data_' Info.Sides{iSide} '.scanData = script_getROIdata(thisView,scanDataVar,glmInfo.analysisBaseNames_Scans,roiNames,glmInfo.analysisScanNum,' q 'overlays' q ');']);
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData = script_getROIdata(thisView,scanDataVar,glmInfo.analysisBaseNames_Scans,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi,glmInfo.analysisScanNum,' q 'overlays' q ');']);
    end
end

%% GLM ROI analysis
% results from GLM roi analysis inform the pRF analysis so perform it first

%% perform ROI analysis
% NOTE: selecting data should happen outside of function
% roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');
for iSide = 1:length(Info.Sides)
    eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
    for iROI = 1:length(ROInames)
        eval(['roidata = data.' Info.Sides{iSide} '.' ROInames{iROI} ';']);
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = script_ROIAnalysis(roidata,Info,glmInfo,stimInfo,plotInfo,subjectInfo,glmInfo.analysisScanNum,' q 'overlays' q ',ROInames{iROI});']);
    end
end


%% save data
% save so don't need to load again
save(saveName,'data','-v7.3');

% first get view so we have the ROIs
thisView = getMLRView;

%% pRF analysis
% perform pRF analysis (restricted to auditory responsive voxels * [3 3 3] sphere (ARexp ROI))
% add analysis per scan
[thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,'ARExp',0);

%% pRF grandient reversals
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase, 'pRF_ARexp','[18 18 21]');
% create ROIs with the names: LeftGR_pRF, RightGR_pRF

%% export pRF overlays
% export group data from volumetric to flatmap space
for iSide = 1:length(subjectInfo.flatmapNames)
    for iGroup = 1:length(glmInfo.groupNames)
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},[pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, '_vol' ],[],[],subjectInfo.flatmapNames{iSide});
        end
    end
end

%% average pRF overlays over depth cortical depth
% pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
overlayNames = cell(size(pRFInfo.analysisNames_Groups));
for iSide = 1:length(subjectInfo.flatmapNames)
    % get  group number and analysis number and get analysis, then get overlays
    % and there names - concat with what makes exported overlay
    for iGroup = 1:length(glmInfo.groupNames)
        groupName = glmInfo.groupNames{iGroup};
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            %             analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
            analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, '_vol' ];
            for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                overlayNames{iGroup}{iAnal}{iOverlay} = [groupName '_' analysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
            end
            thisView = script_averageAcrossDepths(thisView,overlayNames{iGroup}{iAnal},[subjectInfo.flatmapNames{iSide}, 'Volume']);
        end
    end
end

%% get pRF data
% list overlay names and get using script_getOverlayData
overlayNames = cell(size(pRFInfo.analysisNames_Groups));
% overlayNames= ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
% overlayNames = 'averageDepthVol(pRF(PrefCentreFreq,0))';

for iSide = 1:length(subjectInfo.flatmapNames)
    baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
    thisView = viewSet(thisView,'currentbase',baseNum);
    for iGroup = 1:length(glmInfo.groupNames)
        groupName = [];
        groupName = glmInfo.groupNames{iGroup};
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            analysisName = [];
            analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, '_vol' ];
            %             analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
            for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                overlayNames{iGroup}{iAnal}{iOverlay} = ['averageDepthVol(' groupName '_' analysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0))'];
            end
            %                  overlayNames{iGroup}{iAnal}{iOverlay} = [groupName '_' analysisName ' (' pRFOverlayNames{iOverlay} ',0)'];
            
            %             function overlayData = script_getOverlayData(thisView,groupName,analysisName,conditionNames,iScan)
            eval(['data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames{iGroup}{iAnal},[]);'])
            %         for iAnal = 1:length(glmInfo.analysisNames_Groups)
            %             eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,[]);'])
            %         end
            
        end
    end
end

%% Restrict pRF data by ROI
% difference number of analysis in each group
% use the same format as glm - iAnal and loop through for each group

for iSide = 1:length(Info.Sides)
    eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
    for iROI = 1:length(roiNames)
        roi = viewGet(thisView,'roi',roiNames{iROI});
        for iGroup = 1:length(glmInfo.groupNames)
            
            %             for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                %                 analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
                
                analysisName = pRFInfo.analysisNames_Groups{iGroup};
                %                 eval(['dataVar = data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup} ';']);
                eval(['dataVar = data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.overlayData;']);
                
                %         function ROIdata = script_getROIdata(thisView,analysisData,analysisBaseNames,ROI,analysisScanNum,dataType)
                eval(['data.' Info.Sides{iSide} '.' roiNames{iROI} '.' glmInfo.groupNames{iGroup} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} ' = script_getROIdata(thisView,dataVar,analysisName,roi,[],' q 'overlays' q ');']);
            end
        end
    end
end
%% pRF analysis
% create function or add to: voxel comparisions;
% pCF scatter, pCF distribution, pCF correlation
% r = corr2(A,B)

%% get pCF data
pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
for iSide = 1:length(Info.Sides)
    eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
    for iROI = 1:length(roiNames)
        %         roi = viewGet(thisView,'roi',roiNames{iROI});
        for iAnala = 1:length(pRFInfo.analysisNames_Groups{1})
            
            eval(['conA_data = data.' Info.Sides{iSide} '.' roiNames{iROI} '.' glmInfo.groupNames{1} '.' pRFInfo.analysisNames_Groups{1}{iAnala} ';']);
            %             for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            for iAnalb = 1:length(pRFInfo.analysisNames_Groups{2})
                
                
                % data = script_pRFROIAnalysis(conA_data,conB_data,dataNames)
                
                eval(['conB_data = data.' Info.Sides{iSide} '.' roiNames{iROI} '.' glmInfo.groupNames{2} '.' pRFInfo.analysisNames_Groups{2}{iAnalb} ';']);
                
                tmpData{iSide}{iROI}{iAnalb} = script_pRFROIAnalysis(conA_data,conB_data,pRFInfo.pRFOverlayNames);
            end
        end
    end
end



%% Save data
save(saveName,'data','-v7.3');


%% difference map
% go to flatmap groups > take overlay from each group > subtrack them from
% each other > install as new overlay

for iSide = 1:length(Info.Sides)
    overlay = cell(size(pRFInfo.analysisNames_Groups));
    for iGroup = 1:length(glmInfo.groupNames)
        
        thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{iGroup});
        
        % 'overlay'
        %    overlay = viewGet(view,'overlay',[overlayNum],[analysisNum])
        %    overlay = viewGet(view,'overlay',overlayNum,[])
        %    overlay = viewGet(view,'overlay',[],analysisNum)
        %    overlay = viewGet(view,'overlay',[],[])
        %    overlay = viewGet(view,'overlay',overlayNum)
        %    overlay = viewGet(view,'overlay')
        
        %% loop over analysis
        % add this in pRF analysis loop?
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            %         analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
            analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, 'Vol' ];
            
            thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
            overlayNum = viewGet(thisView,'overlayNum','PrefCentreFreq');
            overlay{iGroup}{iAnal} = viewGet(thisView,'overlay',overlayNum);
        end
    end
    for iAnal = 1:length(pRFInfo.analysisNames_Groups{2})
        %     analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
        analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, 'Vol' ];
        
        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
        [ thisView , differenceData ] = script_createDifferenceMaps(thisView,overlay{1}{1},overlay{2}{iAnal});
    end
end

%% Cortical Magnification
% look at Juliens script for CM - cortical distances etc
% define line ROIs for revHa revL RevHp
% get cortical distances
% get curvature?
% get pCF estiamte
% plot

%% plot study information
[ data ] = plot_studyInfo(stimInfo, glmInfo, pRFInfo, Info, plotInfo)

%% quit mrLoadRet
mrQuit()