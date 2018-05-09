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
dataMount = '/Volumes/DataDisk/data';

%% Get subject info
subjectInfo = get_SubjectInfo_CM(iSub);
% Subject ID, flatmap names
saveName = [subjectInfo.subjectID '_data.mat'];

%% move to subject folder, delete any current views, open mrLoadRet and get its view
cd(fullfile(dataMount,Info.studyDir,subjectInfo.subjectID));
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
sHL_organiseData(Info, subjectInfo);
% import, convert and move subject data

%% pre-process
sHL_preprocess(Info, subjectInfo, 0);
% distortion correct
% linear alignment
% non-linear alignment

%% Setup mrLoadRet
% smooth = 1;
[thisView, concatedate] = script_setupmrLoadRet(Info,subjectInfo,glmInfo);
% initiate mrLoadRet
% motion corerection
% smooth (optional)
% group data
% concatenate runs

%% open View
mrLoadRet

%% Import anatomy
thisView = script_importAnatomy(thisView);
% load in:
% reference EPI
% High resolution in-plane T2*
% surfaces
% create flatmaps
% MAKE ORIGINAL FLAT MAPS (LEFT AND RIGHT) USING MAKEFLAT AND NAME THEM
% [freeSurferName{iSubj} '_left_Flat.off'] AND [freeSurferName{iSubj} '_right_Flat.off']
% rotate flatmaps for easy viewing (do before exporting to flatmap space)

%% Tonotopic analysis

%% GLM analysis
% first do Box Car
thisView = script_glmAnalysis(thisView,glmInfo,{'hrfBoxcar'},1);
% HRF = double gamma and box car
% All stims and 8 bins
% don't need to perform weighted mean on individual runs

%% Make flatmaps
% radius = 55
% centre on HG (use R2 and f-test to guide)

%% GLM grandient reversals
% rotate left flatmap 230 and right 290
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,glmInfo.groupNames{1}, glmInfo.analysisNames_Groups{1},'[18 18 21]');

%% ROI CREATION
% create ROIs with the names: 
% LeftAC_glmbc, RightAC_glmbc using gradient reversals - output 4 with alpha overlay output 6  
%   go to each flat base vol, define large ROI around HG, project through 
%   depths (ROIs>transform>expandROI([1 1 6])(post fix _exp to rio name),

%% now create pRF restrict ROI in flat space and project through depths,
% 
% first get view so we have the ROIS
thisView = getMLRView;

% go to each flat base vol, define large ROI around HG on slice 6, project through depths (ROIs>transform>expandROI([30 30 6])(post fix 'ex' to name)),

% eval(['glmbc_rois = [Info.' Info.Sides{iSide} 'ROInames ' q 'ex' q '];']);
% pRFrois = {'LeftpRFrestrict','RightpRFrestrict'}; % this should be defined in setup
glmbc_rois = {'LeftAC_glmbc_vol','RightAC_glmbc_vol'}; % this should be defined in setup

for iSide = 1:length(glmbc_rois)
% roi = viewGet(thisView,'roi','leftpRFrestrict');
% roi = viewGet(thisView,'roi','RightAC');       
baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
thisView = viewSet(thisView,'currentbase',baseNum);
thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide} 'Volume']);

refreshMLRDisplay(thisView.viewNum);
roi = viewGet(thisView,'roi',glmbc_rois{iSide});
outputRoi = convertFromFlatVolumeToBase(roi);
thisView = viewSet(thisView,'newROI',outputRoi);
% post fix Vol to name
end

%% HRF estimate
thisView = script_hrfAnalysis(thisView,glmInfo.groupNames{2});
% estimate hrf using deconvolution

% for iSide = 1:length(glmbc_rois)
%     roiName = [glmbc_rois{iSide} '_vol'];
%     [ x_doubleGamma{iSide}, x_Gamma{iSide}, x_dGamma{iSide} ] = script_hrfROIAnalysis(thisView,roiName,glmInfo);
%     % get data from analysis
%     % use ROI to restrict
%     % perform ROI analysis - average hrf estimate
%     % output result
% end
thisView = getMLRView;
% save result to data
[ data.hrf.x_doubleGamma, data.hrf.x_Gamma, data.hrf.x_dGamma, data.hrf.deconv] = script_hrfROIAnalysis(thisView,'AC_glmbc_vol',glmInfo);
glmInfo.hrfParamsDoubleGamma = data.hrf.x_doubleGamma;
pRFInfo.hrfParamsGamma = data.hrf.x_Gamma;
pRFInfo.hrfParamsDiffofGamma = data.hrf.x_dGamma;

% overlayData = script_getOverlayData(thisView)

%% get av HRF estimate for GLM BOXCAR Gradient Reversals ROI
% use results for GLM
% add option for hrf model
thisView = script_glmAnalysis(thisView,glmInfo,{'hrfDoubleGamma'},1);

%% GLM grandient reversals
% rotate left flatmap 230 and right 290
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo,Info.gradReversalInfo.groupBase, glmInfo.analysisNames_Groups{2},'[18 18 21]');

%% ROI CREATION
% create ROIs with the names: 
% LeftAC_glmdg, RightAC_glmdg using gradient reversals - output 4 with alpha overlay output 6
% LeftRestrict, RightRestrict 
%   go to each flat base vol, define large ROI around HG, project through 
%   depths (ROIs>transform>expandROI([1 1 6])(replace)),

%% Convert GLM data to flatmap space and average over cortical depth
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
%         not sure if I need below
%         baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
%         thisView = viewSet(thisView,'currentbase',baseNum);
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
% for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
%     analysisName = glmInfo.analysisNames_Scans{iAnal};
%     conditionNames{iAnal} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
%     
% end
analysisName = glmInfo.analysisNames_Scans{1};
conditionNames{1} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
analysisName = glmInfo.analysisNames_Scans{3};
conditionNames{2} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);


% save condition names
data.conditions = conditionNames;

% get data from SCANS

analysisName = 'combineTransformOverlays';
for iSide = 1:length(subjectInfo.flatmapNames)
    groupName = [subjectInfo.flatmapNames{iSide} 'Volume'];
    thisView = viewSet(thisView,'curgroup',groupName);
    
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
    % if isempty(overlays)
    %     analysisData = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
    %     overlays = 1:length(analysisData.overlays);
    % end
    % if ~isempty(iScan)
    %     for iCon = 1:length(conditionNames)
    %         overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' analysisName ' (' conditionNames{iCon} ',0))'];
    %     end
    % overlayData = get_overlayData(thisView,overlayNames);
    % else
    for iScan = 1:glmInfo.nScans
        for iAnal = 1:length(glmInfo.hrfModel)*length(glmInfo.nStim)
            overlayNames = [];
            if glmInfo.analysisNStim{iAnal} == length(conditionNames{1});
                for iCon = 1:length(conditionNames{1})
                    overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisBaseNames_Scans{iAnal} '_Scan_' mat2str(iScan) ' (' conditionNames{1}{iCon} ',0))'];
                end
            else
                for iCon = 1:length(conditionNames{2})
                    overlayNames{iCon} = ['averageDepthVol(Scan ' mat2str(iScan) ' - ' glmInfo.analysisBaseNames_Scans{iAnal} '_Scan_' mat2str(iScan)  ' (' conditionNames{2}{iCon} ',0))'];
                end
            end
            overlayData = get_overlayData(thisView,overlayNames);
            eval(['data.' Info.Sides{iSide}, '.scanData.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} =  overlayData']);
            %             eval(['data.' Info.Sides{iSide}, '.scanData.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,iScan);'])
            %             eval(['data.' Info.Sides{iSide}, '.scanData.', glmInfo.analysisBaseNames_Scans{iAnal}, '.overlayData{iScan} = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',conditionNames{iAnal},iScan);'])
        end
    end
end

% get data from GROUPs
for iGroup = 1:length(glmInfo.groupNames)
    %         overlayNames= ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
    for iAnal = 1:length(glmInfo.analysisNames_Groups)
    overlayNames = [];
    conNamesString = [];
    for iCon =1:length(conditionNames{1})
        overlayNames{iCon} = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames_Groups{iAnal} ' (' conditionNames{1}{iCon} ',0))'];
        
        if iCon == 1
            conNamesString  = [conNamesString conditionNames{1}{iCon}];
        else
            conNamesString  = [conNamesString ',' conditionNames{1}{iCon}];
        end
    end
    centriodOverlayName = ['averageDepthVol(' glmInfo.groupNames{iGroup} '_' glmInfo.analysisNames_Groups{iAnal} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
    for iSide = 1:length(subjectInfo.flatmapNames)
        baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
        thisView = viewSet(thisView,'currentbase',baseNum);
        for iAnal = 1:length(glmInfo.analysisNames_Groups)
            eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.overlayData = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',overlayNames,[]);'])
            eval(['data.' Info.Sides{iSide}, '.' glmInfo.groupNames{iGroup} '.', glmInfo.analysisNames_Groups{iAnal}, '.glmCentriod = script_getOverlayData(thisView,[subjectInfo.flatmapNames{iSide},' q 'Volume' q '],' q 'combineTransformOverlays' q ',centriodOverlayName,[]);'])
        end
    end
    end
end

%% Restrict data by ROIs
% Now we need to restrict the data by the ROIs

% first get view so we have the ROIS
thisView = getMLRView;

% loop - save (side), selet:roi (side) and base (side), get data from: side.scans.anal.overlayData.data
% i.e. get data from left flatmap (set: left roi, left group, left base) save this data to left struct

% GROUPS
for iSide = 1:length(Info.Sides)
    
    eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
    baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
    thisView = viewSet(thisView,'currentbase',baseNum);
    
    ROIdata = struct;
    q = char(39);
    for iROI = 1:length(ROInames)
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = struct;']);
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
        
        
        for iGroup = 1:length(glmInfo.groupNames)
            %         eval(['ROI_data_' Info.Sides{iSide} '.' 'glmInfo.groupNames{iGroup} ' = script_getROIdata(thisView,groupDataVar,glmInfo.analysisNames_Groups,roiNames,[],' q 'overlays' q ');']);
            for iAnal = 1:length(glmInfo.analysisNames_Groups)
                
                eval(['groupDataVar = data.' Info.Sides{iSide} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} '.overlayData;']);
                analysisName = glmInfo.analysisNames_Groups{iAnal};
                eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.' glmInfo.groupNames{iGroup} '.' glmInfo.analysisNames_Groups{iAnal} ' = script_getROIdata(thisView,groupDataVar,analysisName,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi,[],' q 'overlays' q ');']);
            end
        end
        
        % SCANS - get data from ROIs
        eval(['scanDataVar = data.' Info.Sides{iSide} '.scanData;']);
        %     eval(['ROI_data_' Info.Sides{iSide} '.scanData = script_getROIdata(thisView,scanDataVar,glmInfo.analysisBaseNames_Scans,roiNames,glmInfo.analysisScanNum,' q 'overlays' q ');']);
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} '.scanData = script_getROIdata(thisView,scanDataVar,glmInfo.analysisBaseNames_Scans,data.' Info.Sides{iSide} '.' ROInames{iROI} '.roi,glmInfo.analysisScanNum,' q 'overlays' q ');']);
    end
end

%% GLM ROI analysis
% results from GLM roi analysis inform the pRF analysis so perform it first
%% perform ROI analysis
% save so don't need to load again
% compare binning for glm to averaging betas
% roiAnalysis = script_ROIAnalysis(roiData,glmInfo.analysisBaseNames_Scans,Info,stimInfo,plotInfo,Info.conditionRunIndex,glmInfo.analysisScanNum,'GLM');
for iSide = 1:length(Info.Sides)
    eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
    for iROI = 1:length(ROInames)
        eval(['roidata = data.' Info.Sides{iSide} '.' ROInames{iROI} ';']);
        
%         e = estimate.hdr;
% [index, threshold] = cal_R2threshold(r2data{:}(volumeIndices));
% e = e(:,:,index(:));
% could add index option to script_ROIfunction or save it to data struct


        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = script_ROIAnalysis(roidata,Info,glmInfo,stimInfo,plotInfo,subjectInfo,glmInfo.analysisScanNum,' q 'overlays' q ',ROInames{iROI});']);
    end
end


%% save data
save(saveName,'data','-v7.3');


%% now create pRF restrict ROI in flat space and project through depths,

% first get view so we have the ROIS
thisView = getMLRView;

% go to each flat base vol, define large ROI around HG on slice 6, project through depths (ROIs>transform>expandROI([30 30 6])(post fix 'ex' to name)),

eval(['pRFrois = [Info.' Info.Sides{iSide} 'ROInames ' q 'ex' q '];']);
% pRFrois = {'LeftpRFrestrict','RightpRFrestrict'}; % this should be defined in setup
pRFInfo.pRFrois = {'LeftAC_glmbc_ex','RightAC_glmbc_ex'}; % this should be defined in setup

for iSide = 1:length(pRFrois)
% roi = viewGet(thisView,'roi','leftpRFrestrict');
% roi = viewGet(thisView,'roi','RightAC');       
baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
thisView = viewSet(thisView,'currentbase',baseNum);
thisView = viewSet(thisView,'curgroup',[subjectInfo.flatmapNames{iSide} 'Volume']);

refreshMLRDisplay(thisView.viewNum);
roi = viewGet(thisView,'roi',pRFrois{iSide});
outputRoi = convertFromFlatVolumeToBase(roi);
thisView = viewSet(thisView,'newROI',outputRoi);
% post fix Vol to name
end
% then convert back to vol space, THEN:  ROIs>combine>(choose both ROIS),action:Union, newName:pRFrestrict
% name = pRFrestrict

%% pRF analysis
% maybe move to before exporting to flatmap and then run glm informed pRF
% later?
% change export and average over depth functiosn to work on indivudal
% overlays as well as all in analysis
for iSide = 1:length(Info.Sides)    
    eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
% get from glm analysis
% average tuning curve sigma
% hrf estimate
[thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,[ pRFInfo.pRFrois{iSide} , '_vol' ],0);
end


%% pRF grandient reversals
% doesn't work?
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);
% create ROIs with the names: pRFleft, pRFright

%% export pRF overlays
% export group data
for iGroup = 1:length(glmInfo.groupNames)
    for iSide = 1:length(subjectInfo.flatmapNames)
        %         baseNum = viewGet(thisView,'baseNum',[subjectInfo.flatmapNames{iSide} 'Volume']);
        %         thisView = viewSet(thisView,'currentbase',baseNum);
%         thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},[pRFInfo.analysisNames_Groups{1}{1}, '_', pRFInfo.pRFrois{iSide}, '_vol' ],[],[],subjectInfo.flatmapNames{iSide});
     
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            %             [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
            % export of flatmap space
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},[pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFInfo.pRFrois{iSide}, '_vol' ],[],[],subjectInfo.flatmapNames{iSide});
        end
    end
end
%% average pRF overlays over depth
% average over cortical depth

% pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
overlayNames = cell(size(pRFInfo.analysisNames_Groups));
for iSide = 1:length(subjectInfo.flatmapNames)
    %     thisView = script_averageAcrossDepths(thisView,overlays,groupName)
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

    %         overlayNames= ['averageDepthVol(' glmInfo.groupNames{iGroup} ' (Ouput 3 - weightedMeanStd(' conNamesString '),0))'];
    
%     overlayNames = 'averageDepthVol(pRF(PrefCentreFreq,0))';

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

%% get weighted mean data
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
%% plot study information
[ data ] = plot_studyInfo(stimInfo, glmInfo, pRFInfo, Info, plotInfo)


%% quit mrLoadRet
mrQuit()