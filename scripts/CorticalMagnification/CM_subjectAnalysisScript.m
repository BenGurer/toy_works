%% CM_subjectAnalysisScript
% Scripted analysis for Comparisions at 7T and Cortical Magnification (CM) studies
% study dates: 2015 - 2018
% Aims:
% Measure cortical magnification
% Compare acquisition protocols
% Compare analysis methods

% load data
% GLM analysis
% pRF analysis
% Gradient reversals and ROI definition
% Cortical distance
% export data from volume to flatmap space
% avearge data over cortical depths
% save data to structure for group analysis

%% save structure
% data.roi.dataset.analysis.voxelParameter
% ie. data.Left.scan1.glm_boxcar_nCons8.R2


%% Get study parameters
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects

% use is pc to set data directory - could do in cm_setupStduyparams
% Info.dataDir
Info.dataDir = '/Volumes/data_PSY/OneDrive - The University of Nottingham/data';
% Info.dataDir = 'E:\data';
q = char(39);

doAntomy = 0;
doPSIR = 0;
doGLMbc = 1;
doMakeFlatmaps = 0;
doHRF = 1;
doGLMdg = 1;
    doConvertOverlays = 1;
        doDeleteOverlays = 1;
        doConvertvol2FlatAvDepth_GLM = 1;
    doGradientReversal_GLM = 1;
    doROIsGLM = 0;
    doGetDATA_GLM = 1;
    doROIRestrict_GLM = 1;
dopRF = 1;
    doConvertvol2FlatAvDepth_pRF = 1;
    doGradientReversal_pRF = 1;
    doROIspRF = 0;
    doGetDATA_pRF = 1;
    doROIRestrict_pRF = 1;
runCorticalMagnification = 1;

%% define subject
iSubs2Run = [1,2,3,4,5,6,7,8];
% iSub= 5;

for iSub = 1:length(iSubs2Run)
    
    %% Get subject info
    subjectInfo = get_SubjectInfo_CM(iSubs2Run(iSub));
    
    % Subject ID, flatmap names
    saveName = [subjectInfo.subjectID '_data.mat'];
    
    %% move to subject folder, delete any current views, open mrLoadRet and get its view
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
    % deleteView(thisView);
    
    %% either load data or mrView
    if exist(saveName, 'file')
        load(saveName)
        if isfield(data, 'hrf')
            glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
            pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
            pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;
        end
    else
        data = struct;
    end
    
    %% open View
    mrLoadRet
    
    thisView = getMLRView;
    
    refreshMLRDisplay(thisView.viewNum);
    
    %% Import anatomy
    if doAntomy
        thisView = script_importAnatomy(thisView,Info,subjectInfo);
        % load in:
        % reference EPI
        % High resolution in-plane T2*
        % surfaces
        
    end
    
    %% Import PSIR
    if doPSIR
        % rotate surfaces and then get view
        % save rotations to subjectInfo function
        thisView = getMLRView;
        
        % Import PSIR
        thisView = importPSIRasOverlay(thisView,Info,subjectInfo);
        
        % import flatmaps without FNIRT coords (after flatmaps have been created)
        thisView = script_importFlatmaps(thisView,Info,subjectInfo);
        
    end
    
    %% Tonotopic analysis
    if doGLMbc
        %% GLM Analysis with Boxcar HRF
        % Perform GLM analysis with boxcar model of HRF
        % hrf = Box Car
        % concatenated runs only
        % [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)
        thisView = script_glmAnalysis(thisView,glmInfo,glmInfo.groupNames,{'hrfBoxcar'},1,1);
        
        %% ROI CREATION
        % create ROIs with the names:
        % LeftAR and RightAR: Group=Sparse; Analysis=glm_hrfboxcar; Overlay=f-test (FDR adjusted)- set min to 0.05
        % Create ROI - continuous voxels;
        % ROIs>transform>expandROI([3 3 3])(convolves ROI with a sphere with a diameter of 3^3 voxels)
        % name=(Left or Right)ARexp
        % Project through depths 0.3 to 0.7 to remove voxels outside of grey matter for ALL ROIs
        % combine LeftAR and RightAR = AR & combine LeftARexp and RightARexp = ARexp
        
        disp('create ROIs with the names:')
        disp('LeftAR and RightAR: Group=Sparse; Analysis=glm_hrfboxcar; Overlay=f-test (FDR adjusted)- set min to 0.05')
        disp('Create ROI - continuous voxels;')
        
        keyboard
        
        %%%%%%!!!!!!!!!%%%%%%%%%%%%
        % add combine to GR rois
            roiList(1) = viewGet(thisView,'roinum','LeftAR');
            
            roiList(2) = viewGet(thisView,'roinum','RightAR');            
            
            % get rois
            rois = viewGet(thisView,'roi',roiList);
            
            % Project ROIs across all cortical depths
            for iSide = 1:3
                clear outputRoi
                outputRoi = expandROI(rois(iSide),[3 3 3],[]);
                outputRoi.name = [outputRoi.name '_ex'];
                thisView = viewSet(thisView,'newROI',outputRoi);
            end
            
            thisView = convertROICorticalDepth(thisView,params,varargin);
            
           thisView = combineROIs(thisView,roi1,roi2,action,newName);
        %%%%%%!!!!!!!!!%%%%%%%%%%%%
    end
    
    %% Make flatmaps
    if doMakeFlatmaps
        keyboard
        % radius = 55
        % centre on HG (use R2 and f-test to guide)
        % use default names
        % resolution = 3; method = mrFlatMesh
        % rotate flatmaps for easy viewing (do before exporting to flatmap space)
    end
    
    
    % get view so we have the ROIs
    thisView = getMLRView;
    
    %% HRF estimate
    % estimate hrf using deconvolution
    if doHRF
        thisView = script_hrfAnalysis(thisView,glmInfo.groupNames{2});
        
        %% get av HRF estimate for Auditory Responsive (AR) ROI
        % use results for GLM
        % save result to data
        [ thisView, data.hrf.x_doubleGamma, data.hrf.x_Gamma, data.hrf.x_dGamma, data.hrf.estimate, data.hrf.deconv, data.hrf.deconvTW] = script_hrfROIAnalysis(thisView,'AR',glmInfo);
        % save to analysis structures to pass to functions
        glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
        pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
        pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;
    end
    % Estimate and fit average ROI tuning width using recentre method from deconvolution data
    % [ thisView, x_g, x_r, tw_Deconv, estimate, ~, nVoxels] = script_centredTWROIAnalysis(thisView,'AR',glmInfo);
    
    %% GLM Analysis with Double Gamma HRF
    if doGLMdg
        % Perform GLM analysis with double gamma model of HRF
        % hrf = Box Car
        % All runs
        % [thisView] = script_glmAnalysis(thisView,glmInfo,groupNames,hrfModel,runSplitHalf,runTonotopic)        
        glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
        thisView = script_glmAnalysis(thisView,glmInfo,glmInfo.groupNames,{'hrfDoubleGamma'},1,1);
    end
    
    %% get condition names
    getConditionNames = cell(1,length(glmInfo.nStim));
    % save condition names
    if exist('data') && isfield(data, 'conditions')
        conditionNames = data.conditions;
    else
        analysisName = glmInfo.analysisNames_Scans{1};
        conditionNames{1} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
        analysisName = glmInfo.analysisNames_Scans{3};
        conditionNames{2} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
        data.conditions = conditionNames;
    end
    
    %% Convert overlays
    % convert pCF overlays to NERB
    % all tonotopic estimates - groups and scans
    if doConvertOverlays
        thisView = getMLRView;
        for iGroup = 1:length(glmInfo.groupNames)
            % select group
            groupName = glmInfo.groupNames{iGroup};
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
                thisView = viewSet(thisView,'curgroup',groupName);
            end
            % Select analysis
            for iAnal = 1:length(glmInfo.analysisNames)
                analysisName = glmInfo.analysisNames{iAnal};
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                % get condition names and create concatenate string
                conNamesString = [];
                for iCon = 1:length(conditionNames{1})
                    	stimOne = stimInfo.stimNERBs(1);
                        stimTwo = stimInfo.stimNERBs(2);
                    if iCon == 1
                        conNamesString  = [conNamesString conditionNames{1}{iCon}];
                    else
                        conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                    end
                end
                
                % glmInfo.voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
                
                overlayNum = zeros(1,length(glmInfo.voxelPropertyNames));
                for i = 1:length(glmInfo.voxelPropertyNames)-1
                    overlayNum(i) = viewGet(thisView,'overlayNum',['Ouput ' num2str(i) ' - weightedMeanStd_CM(' conNamesString ')']);
                end
                overlayNum(end) = viewGet(thisView,'overlayNum',['Ouput 1 - indexMax(' conNamesString ')']);
                
                for i = 1:length(overlayNum)
                    overlayIN = viewGet(thisView,'overlay',overlayNum(i));
                    [ thisView , ~ ] = convertOverlay_GLMCF2NERB(thisView,overlayIN,stimOne,stimTwo,[glmInfo.voxelPropertyNames{i} '_nERB']);
                end
                
            end
            
        end
        %% Scans
        % Select group
        groupName = glmInfo.scanGroupName;
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        for iScan = 1:glmInfo.nScans
            thisView = viewSet(thisView,'curScan', iScan);
            for iAnal = 1:length(glmInfo.analysisNames_nCons)
                % select analysis
                analysisName = [glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan)];
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                % get condition names and create concatenate string
                if glmInfo.analysisNStim{iAnal} == length(conditionNames{1})
                    overlayFlatNames = cell(1,length(conditionNames{1}));
                    conNamesString = [];
                    for iCon = 1:length(conditionNames{1})
                     	stimOne = stimInfo.stimNERBs(1);
                        stimTwo = stimInfo.stimNERBs(2);
                        if iCon == 1
                            conNamesString  = [conNamesString conditionNames{1}{iCon}];
                        else
                            conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                        end
                    end
                else
                    stimOne = stimInfo.stimNERBs_bin(1);
                    stimTwo = stimInfo.stimNERBs_bin(2);
                    conNamesString = [];
                    for iCon = 1:length(conditionNames{2})
                        if iCon == 1
                            conNamesString  = [conNamesString conditionNames{2}{iCon}];
                        else
                            conNamesString  = [conNamesString ',' conditionNames{2}{iCon}];
                        end
                    end
                end
                
                overlayNum = zeros(1,length(glmInfo.voxelPropertyNames));
                for i = 1:length(glmInfo.voxelPropertyNames)-1
                    overlayNum(i) = viewGet(thisView,'overlayNum',['Ouput ' num2str(i) ' - weightedMeanStd_CM(' conNamesString ')']);
                end
                overlayNum(end) = viewGet(thisView,'overlayNum',['Ouput 1 - indexMax(' conNamesString ')']);
                
                for i = 1:length(overlayNum)
                    overlayIN = viewGet(thisView,'overlay',overlayNum(i));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     change input below to be stim nERB 1 and 2 rather than stimInfo
                    [ thisView , ~ ] = convertOverlay_GLMCF2NERB(thisView,overlayIN,stimOne,stimTwo,[glmInfo.voxelPropertyNames{i} '_nERB']);
                end
            end
        end
        mrSaveView(thisView);
    end
    %% Convert GLM data to flatmap space and average over cortical depth
    % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
    % export scan data
    % delete overlays if no longer needed
    if doConvertvol2FlatAvDepth_GLM
        
        % auto delete all overlay sbecause its a pain to do
        if doDeleteOverlays
            thisView = getMLRView;
            for iSide = 1:length(subjectInfo.flatmapNames)
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                overlayNumbers = 1:length(analysisData.overlays);
                thisView = viewSet(thisView,'deleteoverlay',overlayNumbers);
            end
        end
        
        
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
            thisView = script_averageAcrossDepths(thisView,[],[subjectInfo.flatmapNames{iSide}, 'Volume'],1);
        end
    
    % save view
    mrSaveView(thisView);
    
    end
    
    %% GLM grandient reversals
    % calculate gradient reversals using GLM (double gamma) analysis.
    % pCF estimation = Juliensdebiased method
    if doGradientReversal_GLM
        groupName = glmInfo.groupNames{1}; % 'ConcatenationSparse';
        if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
            thisView = viewSet(thisView,'curgroup',groupName);
        end
        analysisName = glmInfo.analysisNames{2}; % 'glm_hrfDoubleGamma';
        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
        conNamesString = [];
        for iCon = 1:length(conditionNames{1})
            if iCon == 1
                conNamesString  = [conNamesString conditionNames{1}{iCon}];
            else
                conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
            end
        end
        
%         overlayNum = viewGet(thisView,'overlayNum',['Ouput 3 - weightedMeanStd_CM(' conNamesString ')']);        
        overlayNum = viewGet(thisView,'overlayNum','julien_pCF_nERB');
        %     overlayIN = viewGet(thisView,'overlay',overlayNum);
        thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase, glmInfo.analysisNames_Groups{2},overlayNum,'[18 18 21]');
       
        mrSaveView(thisView)
    end
    
    %% ROI CREATION
    if doROIsGLM
    disp('create ROIs with the names: LeftGR_GLM, LeftGRa_GLM, LeftGRp_GLM, RightGR_GLM, RightGRa_GLM, RightGRp_GLM based on gradient reversals, unsmoothed tonotopic maps and f-test maps')
    disp('Also, line ROIs for each  reversal with the names: LeftHighRevA_GLM, LeftLowRev_GLM, LeftHighRevP_GLM, RightHighRevA_GLM, RightLowRev_GLM, RightHighRevP_GLM')
    %%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%%%%%%%%%%%%%
    % automate combining ROIS
    
    keyboard
    % create ROIs with the names:
    % LeftGR_GLM, LeftGRa_GLM, LeftGRp_GLM, RightGR_GLM, RightGRa_GLM, RightGRp_GLM based on gradient reversals, unsmoothed tonotopic maps and f-test maps
    
    % Also, line ROIs for each  reversal with the names:
    % LeftHighRevA_GLM, LeftLowRev_GLM, LeftHighRevP_GLM, RightHighRevA_GLM, RightLowRev_GLM, RightHighRevP_GLM
    
    %% get GLM data
    % data has now be converted to flatmap space and averaged across cortical depth,
    % the following gets data from flatmap group, in flatmap space, from one depth (middle) that is the average across cortical depth.
    % create names to get data from overlays and save using structure side.Group.anal.data{iScan}
    end
    if doGetDATA_GLM
        %% get data from SCANS
        analysisName = 'combineTransformOverlays';
        % voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
        for iSide = 1:length(subjectInfo.flatmapNames)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            groupName = [subjectInfo.flatmapNames{iSide} 'Volume'];
            thisView = viewSet(thisView,'curgroup',groupName);
            thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
            for iScan = 1:glmInfo.nScans
                for iAnal = 1:length(glmInfo.analysisNames_nCons)
                    % define overlay names
                    % First, deteremine how many stimuli
                    if glmInfo.analysisNStim{iAnal} == length(conditionNames{1})
                        overlayFlatNames = cell(1,length(conditionNames{1}));
                        conNamesString = [];
                        for iCon = 1:length(conditionNames{1})
                            overlayFlatNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' conditionNames{1}{iCon} ',0))'];
                            if iCon == 1
                                conNamesString  = [conNamesString conditionNames{1}{iCon}];
                            else
                                conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                            end
                        end
                    else
                        overlayFlatNames = cell(1,length(conditionNames{2}));
                        conNamesString = [];
                        for iCon = 1:length(conditionNames{2})
                            overlayFlatNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan)  ' (' conditionNames{2}{iCon} ',0))'];
                            if iCon == 1
                                conNamesString  = [conNamesString conditionNames{2}{iCon}];
                            else
                                conNamesString  = [conNamesString ',' conditionNames{2}{iCon}];
                            end
                        end
                    end
                    
                    % beta weights
                    overlayData = get_overlayData(thisView,overlayFlatNames);
                    %  eval(['data.' Info.Sides{iSide}, '.scans.', glmInfo.analysisNames_nCons{iAnal}, '.betas{iScan} =  overlayData;']);
                    eval(['data.' Info.Sides{iSide}, '.scans_',  num2str(iScan) '.', glmInfo.analysisNames_nCons{iAnal}, '.betas =  overlayData;']);
                    
                    % r2 overlay name
                    r2OverlayName = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (r2,0))'];
                    
                    % voxel property estmate overlay names
                    %% change here if overlay names change
                    voxelPropertyOverlayName = cell(1,length(glmInfo.voxelPropertyNames));
                    %                     for iName = 1:length(glmInfo.voxelPropertyNames) - 1
                    %                         voxelPropertyOverlayName{iName} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (Ouput ' num2str(iName) ' - weightedMeanStd(' conNamesString '),0))'];
                    %                     end
                    %                     voxelPropertyOverlayName{end} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (Ouput 1 - indexMax(' conNamesString '),0))'];
                    %
                    for iName = 1:length(glmInfo.voxelPropertyNames) - 1
                        voxelPropertyOverlayName{iName} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' glmInfo.voxelPropertyNames{iName} '_nERB,0))'];
                    end
                    voxelPropertyOverlayName{end} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisNames_nCons{iAnal} '_Scan_' mat2str(iScan) ' (' glmInfo.voxelPropertyNames{end} '_nERB,0))'];
                    
                    
                    % get overlay data using get_overlayData - outputs a structure with fields: .data & .name
                    % R2
                    clear tempData
                    tempData = get_overlayData(thisView,r2OverlayName);                    
%                     eval(['data.' Info.Sides{iSide}, '.scans.' , glmInfo.analysisNames_nCons{iAnal}, '.r2{iScan}  = tempData;']);
                    eval(['data.' Info.Sides{iSide}, '.scans_', num2str(iScan), '.' , glmInfo.analysisNames_nCons{iAnal}, '.r2  = tempData;']);
                    
                    % voxel property estiamtes
                    for iName = 1:length(glmInfo.voxelPropertyNames)
                        clear tempData
                        tempData = get_overlayData(thisView,voxelPropertyOverlayName{iName});
%                         eval(['data.' Info.Sides{iSide}, '.scans.' , glmInfo.analysisNames_nCons{iAnal}, '.', glmInfo.voxelPropertyNames{iName},'{iScan} = tempData;']);
                        eval(['data.' Info.Sides{iSide}, '.scans_', num2str(iScan), '.' , glmInfo.analysisNames_nCons{iAnal}, '.', glmInfo.voxelPropertyNames{iName},' = tempData;']);
                    end
                end
            end
        end
        
        %% get data from GROUPs
        % glmInfo.voxelPropertyNames = {'Centriod','Spread','julien_pCF','julien_pTW','indexMax'};
        for iSide = 1:length(subjectInfo.flatmapNames)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
            thisView = viewSet(thisView,'currentbase',baseNum);
            for iGroup = 1:length(glmInfo.groupNames)
                for iAnal = 1:length(glmInfo.analysisNames)
                    % define overlay names
                    % betas overlay names
                    overlayFlatNames = cell(1,length(conditionNames{1}));
                    conNamesString = [];
                    for iCon =1:length(conditionNames{1})
                        overlayFlatNames{iCon} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' conditionNames{1}{iCon} ',0))'];
                        if iCon == 1
                            conNamesString  = [conNamesString conditionNames{1}{iCon}];
                        else
                            conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
                        end
                    end
                    
                    % r2 overlay name
                    r2OverlayName = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (r2,0))'];
                    
                    % voxel property estmate overlay names
                    %% change here if overlay names change
%                     voxelPropertyOverlayName = cell(1,length(glmInfo.voxelPropertyNames));
%                     for iName = 1:length(glmInfo.voxelPropertyNames) - 1
%                         voxelPropertyOverlayName{iName} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (Ouput ' num2str(iName) ' - weightedMeanStd(' conNamesString '),0))'];
%                     end
%                     voxelPropertyOverlayName{end} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (Ouput 1 - indexMax(' conNamesString '),0))'];
%                     
                    voxelPropertyOverlayName = cell(1,length(glmInfo.voxelPropertyNames));
                    for iName = 1:length(glmInfo.voxelPropertyNames) - 1
                        voxelPropertyOverlayName{iName} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' glmInfo.voxelPropertyNames{iName} '_nERB,0))'];
                    end
                    voxelPropertyOverlayName{end} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames{iAnal} ' (' glmInfo.voxelPropertyNames{end} '_nERB,0))'];
                    
                    
                    % get overlay data using get_overlayData - outputs a structure with fields: .data & .name
                    % R2
                    clear tempData
                    tempData = get_overlayData(thisView,r2OverlayName);
                    eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.r2 = tempData;']);
                    
                    % beta weights
                    clear tempData
                    tempData = get_overlayData(thisView,overlayFlatNames);
                    eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.betas = tempData;']);
                    
                    % voxel property estiamtes
                    for iName = 1:length(glmInfo.voxelPropertyNames)
                        clear tempData
                        tempData = get_overlayData(thisView,voxelPropertyOverlayName{iName});
                        eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames{iAnal}, '.', glmInfo.voxelPropertyNames{iName},' = tempData;']);
                    end
                    
                end
            end
        end
    end
    
    %% Restrict data by ROIs
    % Now we need to restrict the data by the ROIs
    if doROIRestrict_GLM
        
        % first get view so we have the ROIS
        thisView = getMLRView;
        
        q = char(39);
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
            baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
            thisView = viewSet(thisView,'currentbase',baseNum);
            
            for iROI = 1:length(ROInames)
                %         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = struct;']);
                eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
                
                % Restrict Group data
                for iGroup = 1:length(glmInfo.groupNames)
                    for iAnal = 1:length(glmInfo.analysisNames_Groups)/length(glmInfo.groupNames)
                        analysisName = glmInfo.analysisNames_Groups{iAnal};
                        
%                         % restrict r2
%                         eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.r2.data;']);
%                         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.r2  = get_ROIdata(groupDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
%                         clear groupDataVar
%                         
%                         % restrict betas
%                         eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas.data;']);
%                         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas  = get_ROIdata(groupDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
%                         
%                         clear groupDataVar
%                         
%                         % restrict estimates
%                         for iName = 1:length(voxelPropertyNames)
%                             eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.' voxelPropertyNames{iName} '.data;']);
%                             eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.' voxelPropertyNames{iName} '  = get_ROIdata(groupDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
%                             
%                             clear groupDataVar
%                         end
                        
                        % restrict r2
                        eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.r2.data;']);
                        eval(['data.'  ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.r2  = get_ROIdata(groupDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
                        clear groupDataVar
                        
                        % restrict betas
                        eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas.data;']);
                        eval(['data.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.betas  = get_ROIdata(groupDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
                        
                        clear groupDataVar
                        
                        % restrict estimates
                        for iName = 1:length(glmInfo.voxelPropertyNames)
                            eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} '.data;']);
                            eval(['data.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} '  = get_ROIdata(groupDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
                            
                            clear groupDataVar
                        end
                        
                    end
                end
                
                % Restrict Scan data
                for iAnal = 1:length(glmInfo.analysisBaseNames_Scans)/glmInfo.nScans
                    for iScan = 1:glmInfo.nScans
                        
%                         % restrict r2
%                         tempData = [];
%                         eval(['scanDataVar = data.' Info.Sides{iSide} '.scans.' glmInfo.analysisBaseNames_Scans{iAnal} '.r2{iScan}.data;']);
%                         eval(['tempData = get_ROIdata(scanDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
%                         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData.' glmInfo.analysisBaseNames_Scans{iAnal} '.r2{iScan} = tempData{:};']);
%                         clear scanDataVar
%                         
%                         eval(['scanDataVar = data.' Info.Sides{iSide} '.scans.' glmInfo.analysisBaseNames_Scans{iAnal} '.betas{iScan}.data;']);
%                         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData.' glmInfo.analysisBaseNames_Scans{iAnal} '.betas{iScan} = get_ROIdata(scanDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
%                         clear scanDataVar
%                         
%                         
%                         % restrict estimates
%                         for iName = 1:length(voxelPropertyNames)
%                             
%                             tempData = [];
%                             eval(['scanDataVar = data.' Info.Sides{iSide} '.scans.' glmInfo.analysisBaseNames_Scans{iAnal} '.' voxelPropertyNames{iName} '{iScan}.data;']);
%                             eval(['tempData = get_ROIdata(scanDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
%                             eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData.' glmInfo.analysisBaseNames_Scans{iAnal} '.' voxelPropertyNames{iName} '{iScan} = tempData{:};']);
%                             clear scanDataVar                            
%                         end

                        % restrict r2
                        tempData = [];
                        eval(['scanDataVar = data.' Info.Sides{iSide} '.scans_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.r2.data;']);
                        eval(['tempData = get_ROIdata(scanDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
                        eval(['data.' ROInames{iROI} '.scan_'  num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.r2 = tempData{:};']);
                        clear scanDataVar
                        
                        eval(['scanDataVar = data.' Info.Sides{iSide} '.scans_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.betas.data;']);
                        eval(['data.' ROInames{iROI} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.betas = get_ROIdata(scanDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
                        clear scanDataVar
                        
                        
                        % restrict estimates
                        for iName = 1:length(glmInfo.voxelPropertyNames)
                            
                            tempData = [];
                            eval(['scanDataVar = data.' Info.Sides{iSide} '.scans_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.' glmInfo.voxelPropertyNames{iName} '.data;']);
                            eval(['tempData = get_ROIdata(scanDataVar,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi);']);
                            eval(['data.' ROInames{iROI} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{iAnal} '.' glmInfo.voxelPropertyNames{iName} ' = tempData{:};']);
                            clear scanDataVar                            
                        end
                        
                    end
                end
            end
        end

    
    %% GLM ROI analysis
    % results from GLM roi analysis inform the pRF analysis so perform it first
    
    %% perform ROI analysis
    % NOTE: selecting data should happen outside of function
    % roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');
    % for iSide = 1:length(Info.Sides)
    %     eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
    %     for iROI = 1:length(ROInames)
    %         eval(['roidata = data.' Info.Sides{iSide} '.' ROInames{iROI} ';']);
    %         eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = script_ROIAnalysis(roidata,Info,glmInfo,stimInfo,plotInfo,subjectInfo,glmInfo.analysisScanNum,' q 'overlays' q ',ROInames{iROI});']);
    %     end
    % end
    
    %% check data
    % quick plots to make sure its worked
    for iSide = 1:length(Info.Sides)
        eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
        iROI = 1;
        iGroup = 1;
        iAnal = 2;
        iName = 3;
        
%         eval(['roipCFdataA = data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' glmInfo.analysisNames_Groups{iAnal} '.' voxelPropertyNames{iName} ';']);
%         eval(['roipCFdataB = data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' glmInfo.analysisNames_Groups{iAnal} '.' voxelPropertyNames{iName} ';']);
        eval(['roipCFdataA = data.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} ';']);
        eval(['roipCFdataB = data.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' glmInfo.analysisNames_Groups{iAnal} '.' glmInfo.voxelPropertyNames{iName} ';']);
       
        figure
        subplot(2,1,1)
        histogram(cell2mat(roipCFdataA))
        hold on
        histogram(cell2mat(roipCFdataB))
        
        
%         eval(['roiBetadataA = data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
%         eval(['roiBetadataB = data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
        eval(['roiBetadataA = data.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
        eval(['roiBetadataB = data.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' glmInfo.analysisNames_Groups{iAnal} '.betas;']);
        
        betas_mv_A = cal_movingAverage(cell2mat(roiBetadataA'));
        betas_mv_B = cal_movingAverage(cell2mat(roiBetadataB'));
        
        subplot(2,1,2)
        plot(mean(betas_mv_A,2))
        hold on
        plot(mean(betas_mv_B,2))
        
        for iScan = 1:4
%             eval(['roiSplitBetas{iScan} = data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData.' glmInfo.analysisBaseNames_Scans{4} '.betas{iScan};']);
            eval(['roiSplitBetas{iScan} = data.' ROInames{iROI} '.scan_' num2str(iScan) '.' glmInfo.analysisBaseNames_Scans{4} '.betas;']);
            splitBetas{iScan} = cell2mat(roiSplitBetas{iScan}');
        end
        [splitMeanA, ROI_data, Voxel_data, totalROIpCF] = cal_splitMean(splitBetas{1},splitBetas{3});
        
        [splitMeanB, ROI_data, Voxel_data, totalROIpCF] = cal_splitMean(splitBetas{2},splitBetas{4});
        
        figure
        for i =1:8
            subplot(2,4,i)
            plot(splitMeanA(i,:))
            hold on
            plot(splitMeanB(i,:))
        end
        
    end
        %% save data
    % save so don't need to load again
    save(saveName,'data','-v7.3');
    
    %% save view
    mrSaveView(thisView);
    
    end
    
    if dopRF
        % first get view so we have the ROIs
        thisView = getMLRView;
        
        %% pRF analysis
        % perform pRF analysis (restricted to auditory responsive voxels * [3 3 3] sphere (ARexp ROI))
        % ADD analysis per scan
        % add these to study setup function
        pRFrestrictROI = 'ARexp';
        pRFanalysisName = ['pRF_', pRFrestrictROI];
        pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
        pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;
        
        [thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,pRFrestrictROI,1,0);
        
        %% pRF grandient r     98lk;kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkeversals
        if doGradientReversal_pRF
            thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase,pRFanalysisName,pRFInfo.pRFgradientReversalOverlay,'[18 18 21]');
            % thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,groupBase,analysisBase,overlayNumber,smoothingParams)            
           
        end
        if doROIspRF
            
             keyboard
            % create ROIs with the names:
            % LeftGR_pRF, LeftGRa_pRF, LeftGRp_pRF, RightGR_pRF,
            % RightGRa_pRF, RightGRp_pRF based on gradient reversals, unsmoothed tonotopic maps and f-test maps.
            
            % Also, line ROIs for each  reversal with the names:
            % LeftHighRevA_pRF, LeftLowRev_pRF, LeftHighRevP_pRF, RightHighRevA_pRF, RightLowRev_pRF, RightHighRevP_pRF
            
        end
        
        %% export  pRF overlays and average over depth cortical depth
        % pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
        
        if doConvertvol2FlatAvDepth_pRF
            
            %% Groups
            for iSide = 1:length(subjectInfo.flatmapNames)
                
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                    thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                end
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                end
                
                for iGroup = 1:length(glmInfo.groupNames)
                    groupName = glmInfo.groupNames{iGroup};
                    
                    for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                        pRFanalysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_',  pRFrestrictROI];
                        % get overlay names and numbers
                        overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                        overlayNum = zeros(1,length(pRFInfo.pRFOverlayNames));
                        overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                        for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                            overlayFlatNames{iOverlay} = [groupName '_' pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                            overlayNum(iOverlay) = iOverlay;
                            overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                        end
                        
                        % export group data from volumetric to flatmap space
                        % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
                        thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},pRFanalysisName,[],overlayNum,subjectInfo.flatmapNames{iSide});
                        
                        % average over cortical depth
                        thisView = script_averageAcrossDepths(thisView,overlayFlatNames,[subjectInfo.flatmapNames{iSide}, 'Volume'],1);
                        
                    end
                end
            end
        end
        
        if doGetDATA_pRF
            
            
            %% get group data and restrict by ROIs
            for iSide = 1:length(Info.Sides)
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                    thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                end
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                end
                eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
                thisView = viewSet(thisView,'currentbase',baseNum);
                
                for iGroup = 1:length(glmInfo.groupNames)
                    groupName = glmInfo.groupNames{iGroup};
                    for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                        pRFanalysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_',  pRFrestrictROI];
                        overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                        overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                        for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
%                             % get overlay names
%                             overlayFlatNames{iOverlay} = [groupName '_' pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
%                             overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
%                             clear tempData
%                             tempData = get_overlayData(thisView,overlay2Get{iOverlay});
%                             % get data
%                             eval(['data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, ' = tempData;']);
%                             
%                             % restrict by ROI
%                             for iROI = 1:length(roiNames)
%                                 eval(['roi = data.' Info.Sides{iSide} '.' roiNames{iROI} '.roi;']);
%                                 eval(['data.', Info.Sides{iSide}, '.', roiNames{iROI}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '  = get_ROIdata(tempData.data,roi);']);
%                                 
%                             end

                          % get overlay names
                            overlayFlatNames{iOverlay} = [groupName '_' pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                            overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                            clear tempData
                            tempData = get_overlayData(thisView,overlay2Get{iOverlay});
                            % get data
                            eval(['data.' Info.Sides{iSide}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, ' = tempData;']);
                            
                            % restrict by ROI
                            for iROI = 1:length(roiNames)
                                eval(['roi = data.' Info.Sides{iSide} '.' roiNames{iROI} '.roi;']);
                                eval(['data.', roiNames{iROI}, '.', glmInfo.groupNames{iGroup}, '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '  = get_ROIdata(tempData.data,roi);']);
                                
                            end
                            
                        end
                    end
                end
            end
            
            %% Scans
            for iScan = 1:glmInfo.nScans
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                    thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                end
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                end
                pRFanalysisName = ['pRF_',  pRFrestrictROI, '_Scan - ' num2str(iScan)];
                overlayNames = cell(1,length(pRFInfo.pRFOverlayNames));
                overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                overlayNum = zeros(1,length(pRFInfo.pRFOverlayNames));
                overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                    
                    overlayNames{iOverlay} = [pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                    overlayFlatNames{iOverlay} = ['Scan ' num2str(iScan) ' - '  pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                    overlayNum(iOverlay) = iOverlay;
                    overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                end
                
                % export group data from volumetric to flatmap space
                % [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
                % get overlay number
                for iSide = 1:length(subjectInfo.flatmapNames)
                    thisView = script_covertData2FlatmapSpace(thisView,glmInfo.scanGroupName,pRFanalysisName,iScan,[],subjectInfo.flatmapNames{iSide});
                    
                    % average over cortical depth
                    thisView = script_averageAcrossDepths(thisView,overlayFlatNames,[subjectInfo.flatmapNames{iSide}, 'Volume'],1);
                end
            end
        end
        
        %% get scan data and restrict by ROIs
        if doROIRestrict_pRF
            for iSide = 1:length(Info.Sides)
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                    thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
                end
                if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
                end
                eval(['roiNames = Info.' Info.Sides{iSide} 'ROInames;']);
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
                thisView = viewSet(thisView,'currentbase',baseNum);
                
                for iScan = 1:glmInfo.nScans
                    for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
                        pRFanalysisName = ['pRF_',  pRFrestrictROI, '_Scan - ' num2str(iScan)];
                        overlayFlatNames = cell(1,length(pRFInfo.pRFOverlayNames));
                        overlay2Get = cell(1,length(pRFInfo.pRFOverlayNames));
                        for iOverlay = 1:length(pRFInfo.pRFOverlayNames)
                            % get overlay names
                            overlayFlatNames{iOverlay} = ['Scan ' num2str(iScan) ' - '  pRFanalysisName ' (' pRFInfo.pRFOverlayNames{iOverlay} ',0)'];
                            overlay2Get{iOverlay} = ['averageDepthVol(' overlayFlatNames{iOverlay} ')'];
                            clear tempData
                            tempData = get_overlayData(thisView,overlay2Get{iOverlay});
                            % get data
%                             eval(['data.', Info.Sides{iSide}, '.scanData.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} '{iScan} = tempData;']);
%                             
%                             % restrict by ROI
%                             for iROI = 1:length(roiNames)
%                                 clear tempROIdata
%                                 eval(['roi = data.' Info.Sides{iSide} '.' roiNames{iROI} '.roi;']);
%                                 %                             eval(['data.', Info.Sides{iSide}, '.', roiNames{iROI}, '.scanData.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '{iScan}  = get_ROIdata(tempData.data,roi);']);
%                                 eval(['tempROIdata = get_ROIdata(tempData.data,roi);']);
%                                 
%                                 eval(['data.', Info.Sides{iSide}, '.', roiNames{iROI}, '.scanData.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '{iScan}  = tempROIdata{:};'])
%                             end
                            
                            eval(['data.', Info.Sides{iSide}, '.scan_', num2str(iScan), '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, ' = tempData;']);
                            
                            % restrict by ROI
                            for iROI = 1:length(roiNames)
                                clear tempROIdata
                                eval(['roi = data.' Info.Sides{iSide} '.' roiNames{iROI} '.roi;']);
                                %  eval(['data.', Info.Sides{iSide}, '.', roiNames{iROI}, '.scanData.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '{iScan}  = get_ROIdata(tempData.data,roi);']);
                                eval(['tempROIdata = get_ROIdata(tempData.data,roi);']);
                                
                                eval(['data.', roiNames{iROI}, '.scan_', num2str(iScan), '.', pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '.', pRFInfo.pRFOverlayNames{iOverlay}, '  = tempROIdata{:};'])
                            end
                            
                        end
                    end
                end
            end
        end
        
        %% pRF analysis
        % create function or add to: voxel comparisions;
        % pCF scatter, pCF distribution, pCF correlation
        % r = corr2(A,B)
        
        %% check data
        % quick plots to make sure its worked
        for iSide = 1:length(Info.Sides)
            eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
            iROI = 1;
            iGroup = 1;
            iAnal = 1;
            iOverlay = 2;
            
            pRFanalysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_',  pRFrestrictROI];
            
%             eval(['roipCFdataA = data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
%             eval(['roipCFdataB = data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
            eval(['roipCFdataA = data.' ROInames{iROI} '.' glmInfo.groupNames{1} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
            eval(['roipCFdataB = data.' ROInames{iROI} '.' glmInfo.groupNames{2} '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
         
            figure
            histogram(cell2mat(roipCFdataA))
            hold on
            histogram(cell2mat(roipCFdataB))
            
            for iScan = 1:4
                %                 eval(['roiSplitpCF{iScan} = data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} '{iScan};']);
                eval(['roiSplitpCF{iScan} = data.' ROInames{iROI} '.scan_' num2str(iScan) '.' pRFInfo.analysisNames_Groups{iGroup}{iAnal} '.' pRFInfo.pRFOverlayNames{iOverlay} ';']);
            end
            figure
            subplot(2,1,1)
            scatter(roiSplitpCF{1},roiSplitpCF{3})
            subplot(2,1,2)
            scatter(roiSplitpCF{2},roiSplitpCF{4})
            
            % need to remove nans before - also findout how many
            corrcoef(double(roiSplitpCF{1}'),double(roiSplitpCF{3}'))
            corr(roiSplitpCF{1}',roiSplitpCF{3}')
        end
        
        
    %% save data
    save(saveName,'data','-v7.3');
    
    %% save view
    mrSaveView(thisView);
        
    end
    
    
    if runCorticalMagnification
        %% Cortical Magnification
        
        % get cortical distances
        % get curvature?
        % get pCF estiamte
        % plot
        
        % Define gradient reversal ROIs
        % 	High frequency - line
        % 	Low frequency - line
        % 	Gradient - polygon
        
        % loop over:
        %   sides
        %   rois / anterior/ posterior
        %   group
        %   anlysis
        
        % ROIs could have already been create, if not:
        % Create ROIS
        % keyboard
        % Create ROIs for each gradient reversal and name:
        %   'LeftHa','LeftHp','LeftLow,'LeftGRa''LeftGRp'
        %   'RightHa','RightHp','RightLow','RightGRa''RightGRp'
        %   Post fix with _GLM or _pRF
        
        % run script from here
        
        thisView = getMLRView;
        
        % convert rois from flatmap to volume space
        pRFrestrictROI = 'ARexp';
        pRFanalysisName = ['pRF_', pRFrestrictROI];
        analysisNames = {'glm_hrfDoubleGamma',pRFanalysisName};
        AP = {'a','p'};
        analName = {'GLM', 'pRF'};
        for iSide = 1:length(Info.Sides)
            if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',[subjectInfo.flatmapNames{iSide}, 'Volume'])
                thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide}, 'Volume']);
            end
            if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum','combineTransformOverlays')
                thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum','combineTransformOverlays'));
            end
            for iAnal = 1:length(analysisNames)
                for iAP = 1:length(AP)
                    clear rois
                    
                    % make ROI names and get numbers then put in a list
                    HighROIname = [Info.Sides{iSide}, 'High' AP{iAP} '_' analName{iAnal}];
                    roiList(1) = viewGet(thisView,'roinum',HighROIname);
                    
                    LowROIname = [Info.Sides{iSide}, 'Low', '_' analName{iAnal}];
                    roiList(2) = viewGet(thisView,'roinum',LowROIname);
                    
                    GRROIname = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
                    roiList(3) = viewGet(thisView,'roinum',GRROIname);
                    
                    % get rois
                    rois = viewGet(thisView,'roi',roiList);
                    
                    % Project ROIs across all cortical depths
                    for iROI = 1:3
                        clear outputRoi
                        outputRoi = expandROI(rois(iROI),[1 1 6],[]);
                        outputRoi.name = [outputRoi.name '_ex'];
                        thisView = viewSet(thisView,'newROI',outputRoi);
                    end
                    
                    clear roiList
                    % make ROI names and get numbers then put in a list
                    
                    HighROIname = [Info.Sides{iSide}, 'High' AP{iAP} '_' analName{iAnal} '_ex'];
                    roiList(1) = viewGet(thisView,'roinum',HighROIname);
                    
                    LowROIname = [Info.Sides{iSide}, 'Low', '_' analName{iAnal} '_ex'];
                    roiList(2) = viewGet(thisView,'roinum',LowROIname);
                    
                    GRROIname = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal} '_ex'];
                    roiList(3) = viewGet(thisView,'roinum',GRROIname);
                    
                    clear rois_ex
                    % get rois again
                    rois_ex = viewGet(thisView,'roi',roiList);
                    
                    % convert
                    for iROI = 1:3
                        clear outputRoi
                        outputRoi = convertFromFlatVolumeToBase(rois_ex(iROI));
                        outputRoi.name = [outputRoi.name '_Vol'];
                        thisView = viewSet(thisView,'newROI',outputRoi);
                    end
                    
                    % delete expanded rois we don't need
                    thisView = viewSet(thisView,'deleteroi',roiList);
                end
            end
        end
        
        % set view to FNIRTed flatmaps
        % change view to the group / analysis / overlay of interest **** loop over this ****
        % get undistorted flatmap names and send to analysis
        % Perform AverageDepthVol on overlay of interest - calculate in flat space but export to base space
        % Now run cal_tonotopicMagnification to get data
        glmoverlayNames = {'r2', [glmInfo.voxelPropertyNames{3} '_nERB'], [glmInfo.voxelPropertyNames{4} '_nERB']};
        pRFoverlayNames = {pRFInfo.pRFOverlayNames{1}, pRFInfo.pRFOverlayNames{2}, pRFInfo.pRFOverlayNames{3}};
        estimateNames = {'r2', 'pCF', 'pTW'};
        
        for iSide = 1:length(Info.Sides)
            % non- fnirted flatmap names
            flatmapnames.gm = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_GM.off'];
            flatmapnames.wm = [subjectInfo.freeSurferName '_' Info.sides{iSide} '_WM.off'];
            
            % set base
            
            if viewGet(thisView,'currentbase') ~= viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}])
                baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide}]);
                thisView = viewSet(thisView,'currentbase',baseNum);
            end
            
            for iGroup = 1:length(glmInfo.groupNames) % change change this if we want to calculate for both groups or not
                % select group
                groupName = glmInfo.groupNames{iGroup};
                if viewGet(thisView,'curgroup') ~= viewGet(thisView,'groupNum',groupName)
                    thisView = viewSet(thisView,'curgroup',groupName);
                end
                % Select analysis
                for iAnal = 1:length(analysisNames)
                    analysisName = analysisNames{iAnal};
                    
                    if viewGet(thisView,'curAnalysis') ~= viewGet(thisView,'analysisNum',analysisName)
                        thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
                    end
                    
                    for iOverlay = 1:3
                        clear overlayNum
                        switch analysisName
                            case 'glm_hrfDoubleGamma'
                                overlayNum = viewGet(thisView,'overlayNum',glmoverlayNames{iOverlay});
%                                 overlayNames = glmoverlayNames;
                            case pRFanalysisName
                                overlayNum = viewGet(thisView,'overlayNum',pRFoverlayNames{iOverlay});
%                                 overlayNames = pRFoverlayNames;
                        end
                        
                        % check if overlay has already been averaged
                        clear overlayCheck
                        clear overlayname
                        clear tempOverlayData tempOverlayRoiData overlay2getName
                        switch analysisName
                            case 'glm_hrfDoubleGamma'
                                overlay2check = ['averageDepthVol_' Info.Sides{iSide} '(', glmoverlayNames{iOverlay} ')'];
                                
                            case pRFanalysisName
                                overlay2check = ['averageDepthVol_' Info.Sides{iSide} '(', pRFoverlayNames{iOverlay}, ')'];
                        end
                        
                        overlayCheck = viewGet(thisView,'overlayNum',overlay2check);
                        
                        if isempty(overlayCheck)
                            % average overlays over cortical depths
                            [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlayNum)]);
                            params.combineFunction='averageDepthVol';
                            params.combinationMode = 'Apply function to each overlay';
                            params.baseSpace = 1; % perfrom analysis in base space (flatmap)
                            params.exportToNewGroup = 0; % don't export to volume in a new group
                            params.outputName = ['averageDepthVol_', Info.Sides{iSide}];
                            [thisView,params] = combineTransformOverlays(thisView,params);
                        end
                    end
                    switch analysisName
                        case 'glm_hrfDoubleGamma'
                            pCFoverlayname = ['averageDepthVol_' Info.Sides{iSide} '(', glmoverlayNames{2} ')'];
                            
                        case pRFanalysisName
                            pCFoverlayname = ['averageDepthVol_' Info.Sides{iSide} '(', pRFoverlayNames{2}, ')'];
                    end
                    
                    AvOverlayNum = viewGet(thisView,'overlayNum',pCFoverlayname);
                    thisView = viewSet(thisView,'curOverlay',AvOverlayNum);
                    
                    for iAP = 1:length(AP)
                        % make ROI names and get numbers then put in a list
                        clear roiList
                        HighROIname = [Info.Sides{iSide}, 'High' AP{iAP} '_' analName{iAnal} '_ex_Vol'];
                        roiList(1) = viewGet(thisView,'roinum',HighROIname);
                        
                        LowROIname = [Info.Sides{iSide}, 'Low', '_' analName{iAnal} '_ex_Vol'];
                        roiList(2) = viewGet(thisView,'roinum',LowROIname);
                        
                        GRROIname = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal} '_ex_Vol'];
                        roiList(3) = viewGet(thisView,'roinum',GRROIname);
                        
                        % calculate Tonotopic Magnification
                        clear relativeDistances overlayRoiData pathDistances
                        [relativeDistances, overlayRoiData, pathDistances, verticesScanCoords] = cal_tonotopicMagnification(thisView,roiList,flatmapnames,[]);
                        
                        % save to structure
                        roiSaveName = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
                        eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.relativeDistances = relativeDistances;']);
                        eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pCFcheck = overlayRoiData;']);
                        eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pathDistances = pathDistances;']);     
                        
                        
                        interpMethod = 'linear';
                        for iOverlay = 1:3
                            clear tempOverlayData tempOverlayRoiData overlay2getName
                            switch analysisName
                                case 'glm_hrfDoubleGamma'
                                    overlay2getName = ['averageDepthVol_' Info.Sides{iSide} '(', glmoverlayNames{iOverlay} ')'];
                                    
                                case pRFanalysisName
                                    overlay2getName = ['averageDepthVol_' Info.Sides{iSide} '(', pRFoverlayNames{iOverlay}, ')'];
                            end
                            overlay2get = viewGet(thisView,'overlayNum',overlay2getName);
                            tempOverlayData  = viewGet(thisView,'overlaydata',1,overlay2get);
                            tempOverlayRoiData = interpn(tempOverlayData,verticesScanCoords{3}(1,:),verticesScanCoords{3}(2,:),verticesScanCoords{3}(3,:),interpMethod);
                            eval(['data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.' estimateNames{iOverlay} '= tempOverlayRoiData;']);
                            
                        end
                        
                        
                        
                    end
                end
            end
        end
        
    %% save data
    save(saveName,'data','-v7.3');
    
    %% save view
    mrSaveView(thisView);
    
    end

    

    %% quit mrLoadRet
    mrQuit()
end

%% Tidy data
% cortical magnificaiton


pRFrestrictROI = 'ARexp';
pRFanalysisName = ['pRF_', pRFrestrictROI];
analysisNames = {'glm_hrfDoubleGamma',pRFanalysisName};
analysisSaveName = {'GLM','pRF'};
AP = {'a','p'};
analName = {'GLM', 'pRF'};
CorticalDistance = [];
Frequency = [];
Analysis = [];
ROI = [];
r2 = [];
TuningWidth = [];
for iSide = 1:length(Info.Sides)
    for iGroup = 1:length(glmInfo.groupNames)
        groupName = glmInfo.groupNames{iGroup};
        for iAnal = 1:length(analysisNames)
            
            analysisName = analysisNames{iAnal};
            for iAP = 1:length(AP)
                
                roiSaveName = [Info.Sides{iSide}, 'GR' AP{iAP} '_' analName{iAnal}];
                roiName = [Info.Sides{iSide}, 'GR' AP{iAP}];
                
                eval(['tempCorticalDistance = data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.relativeDistances(2,:);']);
                eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pCF;']);
                eval(['tempFrequencycheck = data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pCFcheck{3};']);
                eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.pTW;']);
                eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName '.tonotopicMagnificaion.r2;']);
                
                %             CorticalDistance =
                %             Frequency =
                
%                 if tempFrequency == tempFrequencycheck
                nVoxels = length(tempFrequency);
                tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
                tempROI = repmat(roiName,nVoxels,1);
                
                r2 = [r2; tempR2'];
                TuningWidth = [TuningWidth; tempTuningWidth'];
                CorticalDistance = [CorticalDistance; tempCorticalDistance'];
                Frequency = [Frequency; tempFrequency'];
                if isempty(Analysis)
                        Analysis = tempAnalysis;
                ROI = tempROI;
                else
                Analysis = char(Analysis,tempAnalysis);
                ROI = char(ROI,tempROI);
                end

                
            end
            
        end
    end
end
% Analysis = Analysis(2:end,:);
% ROI = ROI(2:end,:);

T = table(CorticalDistance,Frequency,...
    r2, TuningWidth,...
    Analysis,ROI,...
    'VariableNames',{'CorticalDistance' 'Frequency' 'r2' 'TuningWidth' 'Analysis' 'ROI'});

writetable(T, [saveName, '_CM.csv'])

%% Comparisions
% what to compare?
% what extra things do I need from subjects

%% Comparisions: analysis
% convert units to match - convert overlays so figures are compareable
% convert stimulus to ERB for glm - pRF already convert - steal that code

%% Comparisions: aquistion
% use best analysis


%% difference map
% use averaged over depth overlays just created to make difference maps

% make difference between groups maps - Sparse vs Continuous
% make difference between analysis maps - GLM vs pRF
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


%% Save/export data for group average
% all overlays = averaged across depths
% pCF
% pTW
% Difference maps
% PSIR - need to do the extra regression stuff - thickness and curvature
% Gradient reversals - will need to perform analysis again and export to volume space

%load MNI single subject flat sampling subject space
ssMNI = {'Colin27','Colin27_flipX'};
flatName = {'_34_122_117_Rad60','_141_153_99';'_138_126_114','_34_122_117'};
mniRotation = [230,330;60,230];
flatWarp = {'',['_invFNIRT_' subjects{iSubj}]};
params.anatFileName = fullfile(dataDir,'Anatomy/freesurfer/subjects/', freeSurferName{iSubj}, 'surfRelax',[freeSurferName{iSubj} '_mprage_pp.nii']);
params.flatRes=3;
for iMNI =1:2
  params.path = fullfile(dataDir,'Anatomy/freesurfer/subjects/',ssMNI{iMNI},'surfRelax');
  for iWarp=1:2
    for iSide=1:2
        params.flatFileName = [ssMNI{iMNI} '_' sides{iSide} '_Flat' flatName{iMNI,iSide} '.off'];
        params.outerCoordsFileName = [ssMNI{iMNI} '_' sides{iSide} '_GM_' freeSurferName{iSubj} flatWarp{iWarp} '.off'];
        params.innerCoordsFileName = [ssMNI{iMNI} '_' sides{iSide} '_WM_' freeSurferName{iSubj} flatWarp{iWarp} '.off'];
        params.curvFileName = [ssMNI{iMNI} '_' sides{iSide} '_Curv.vff'];
        base = importFlatOFF(params);
        base.name = [ssMNI{iMNI} '_' sides{iSide} '_Flat' flatWarp{iWarp}];
        thisView = viewSet(thisView, 'newbase', base);
        thisView = viewSet(thisView,'rotate',mniRotation(iMNI,iSide));
    end
  end
end

%%% anatomy
thisView = viewSet(thisView,'curgroup',psirGroup);
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',psirAnalysis));
thisView = viewSet(thisView,'curOverlay',1);
for iSide = 1:2
  %export data to Colin27
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_' sides{iSide} '_Flat']));
  mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'PSIR.nii']));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_flipX_' sides{iSide} '_Flat']));
  mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'PSIRflipX.nii']));
end

%%% Tuning
%set group
thisView = viewSet(thisView,'curgroup',concatenationGroup);
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',functionalAnalysis));
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1','overlayList=[2 3 4 5 6 7 8]');
params.combineFunction='searchlightTuningWidth';
params.additionalArgs = '[3 3 7],1';
params.baseSpace=1;
params.nOutputOverlays=7;
for iSide = 1:2
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',[freeSurferName{iSubj} '_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj} '_lowres']));
  params.outputName=['searchlightTuningWidth ' sides{iSide} ' flat lowres'];
  [thisView,params] = combineTransformOverlays(thisView,params);
  thisView = viewSet(thisView,'curOverlay',mainOverlays(iSubj,[1:10 10+iSide]));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj}]));
  mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'Tonotopy.nii']));
  thisView = viewSet(thisView,'curbase',viewGet(thisView,'basenum',['Colin27_flipX_' sides{iSide} '_Flat_invFNIRT_' subjects{iSubj}]));
  mrExport2SR(thisView.viewNum,fullfile(dataDir,studyDir,'flatExport',[subjects{iSubj} Sides{iSide} 'TonotopyFlipX.nii']));
end


%% plot study information
[ data ] = plot_studyInfo(stimInfo, glmInfo, pRFInfo, Info, plotInfo);


%% Comparisions
% what to compare?
% what extra things do I need from subjects

%% Comparisions: analysis
% convert units to match - convert overlays so figures are compareable
% convert stimulus to ERB for glm - pRF already convert - steal that code

%% Comparisions: aquistion
% use best analysis

%% Group analysis function will:
% import data
% tidy data
% statiscal analysis: average
% plot