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

%% change plotting - create figures outside plotting functions - set paper size etc and use subplotting
%% make data gramm friendly

%% thoughts/notes
% how to remove outliers
% plot glm vs beta weight as a function of frequency as a group average -
% how does it comepare to SL(normalised) normalise both to remove units

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

saveName = [subjectInfo.subjectID '_data.mat'];

%% either load data or mrView

%% move to subject folder, delete any current views, open mrLoadRet and get its view
cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
% deleteView(thisView);

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

%% GLM analysis
thisView = script_glmAnalysis(thisView,glmInfo);
% HRF = double gamma and box car
% All stims and 8 bins
% don't need to perform weighted mean on individual runs

%% GLM grandient reversals
% rotate left flatmap 230 and right 290
thisView = script_flatMapAnalysis(thisView,Info,subjectInfo);

%% ROI CREATION
% create ROIs with the names: 
% LeftAC, RightAC using gradient reversals - output 4 with alpha overlay output 6
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
for iAnal = 1:length(glmInfo.nStim)*length(glmInfo.hrfModel)
    analysisName = glmInfo.analysisNames_Scans{iAnal};
    conditionNames{iAnal} = get_analysisConditionNames(thisView,analysisName,glmInfo.scanGroupName,1);
    
end

% save condition names
data.conditions = conditionNames;

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
        eval(['data.' Info.Sides{iSide} '.' ROInames{iROI} ' = script_ROIAnalysis(roidata,Info,glmInfo,stimInfo,plotInfo,subjectInfo,glmInfo.analysisScanNum,' q 'overlays' q ',ROInames{iROI});']);
    end
end


%% save data
save(saveName,'data','-v7.3');


%% perfrom GLM overlay analysis
fit = cell(1,length(Info.Sides));
weightingType = {'SL','BOLD'};
thisView = viewSet(thisView,'curgroup',glmInfo.groupNames{2});
thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',glmInfo.analysisNames_Groups{1}));
for iWeight = 1:length(weightingType)
    
    switch weightingType{iWeight}
        
        case 'SL'
            namePrefix =  [weightingType{iWeight} 'weighted '];            
            thisView = script_glmOverlayAnalysis(thisView,2:33,namePrefix,0,data.conditions{1});
        case 'BOLD'
            
            for iSide = 1:length(Info.Sides)
                
                eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
                for iROI = 1:length(ROInames)
                    eval(['fit{iSide} = data.' Info.Sides{iSide} '.' ROInames{iROI} '.concatData.' glmInfo.analysisNames_Groups{1} '.roiAnalysis.fit;']);
                end
                
                namePrefix =  [Info.Sides{iSide} '_' weightingType{iWeight} 'weighted '];
                thisView = script_glmOverlayAnalysis(thisView,2:33,namePrefix,fit{iSide},data.conditions{1});
            end
    end
    
end

%% now create pRF restrict ROI in flat space and project through depths,

% first get view so we have the ROIS
thisView = getMLRView;

% go to each flat base vol, define large ROI around HG on slice 6, project through depths (ROIs>transform>expandROI([30 30 6])(post fix 'ex' to name)),

eval(['pRFrois = [Info.' Info.Sides{iSide} 'ROInames ' q 'ex' q '];']);
% pRFrois = {'LeftpRFrestrict','RightpRFrestrict'}; % this should be defined in setup
pRFrois = {'LeftACex','RightACex'}; % this should be defined in setup

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
% for ipRFroi = 1:length(pRFrois)
fit = cell(1,length(Info.Sides));
for iSide = 1:length(Info.Sides)    
    eval(['ROInames = Info.' Info.Sides{iSide} 'ROInames;']);
    for iROI = 1:length(ROInames)    
%     eval(['fit{iSide} = data.' Info.Sides{iSide} '.' ROInames{iROI} '.concatData.' glmInfo.analysisNames_Groups{1} '.roiAnalysis.fit;']);
     eval(['fit{iSide} = data.' Info.Sides{iSide} '.' ROInames{iROI} '.splitData.' glmInfo.analysisBaseNames_Scans{1} '.roiAnalysis.fit;']);
    end
% get info from glm analysis: BOLD ratio, roi av TW
% get info from glm to inform pRF
% BOLD change between conditions
% average tuning curve sigma
glmInfo.m = fit{iSide}(1);
glmInfo.b = fit{iSide}(2);
[thisView, pRFParams] = script_pRFAnalysis(thisView,pRFInfo,glmInfo,[ pRFrois{iSide} , 'Vol' ],1);
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
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
            %             [thisView, analysisData] = script_covertData2FlatmapSpace(thisView,groupName,analysisName,iScan,overlays,flatmapName)
            % export of flatmap space
            thisView = script_covertData2FlatmapSpace(thisView,glmInfo.groupNames{iGroup},[pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFrois{iSide}, 'Vol' ],[],[],subjectInfo.flatmapNames{iSide});
        end
    end
end
%% average pRF overlays over depth
% average over cortical depth

pRFOverlayNames = {'r2','PrefCentreFreq','rfHalfWidth'};
overlayNames = cell(size(pRFInfo.analysisNames_Groups));
for iSide = 1:length(subjectInfo.flatmapNames)
    %     thisView = script_averageAcrossDepths(thisView,overlays,groupName)
        % get  group number and analysis number and get analysis, then get overlays
    % and there names - concat with what makes exported overlay
    for iGroup = 1:length(glmInfo.groupNames)
        groupName = glmInfo.groupNames{iGroup};
        for iAnal = 1:length(pRFInfo.analysisNames_Groups{iGroup})
%             analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};

            analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFrois{iSide}, 'Vol' ];
            for iOverlay = 1:length(pRFOverlayNames)
            overlayNames{iGroup}{iAnal}{iOverlay} = [groupName '_' analysisName ' (' pRFOverlayNames{iOverlay} ',0)'];
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
            analysisName = [pRFInfo.analysisNames_Groups{iGroup}{iAnal}, '_', pRFrois{iSide}, 'Vol' ];
            %             analysisName = pRFInfo.analysisNames_Groups{iGroup}{iAnal};
            for iOverlay = 1:length(pRFOverlayNames)
                overlayNames{iGroup}{iAnal}{iOverlay} = ['averageDepthVol(' groupName '_' analysisName ' (' pRFOverlayNames{iOverlay} ',0))'];
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
                
                tmpData{iSide}{iROI}{iAnalb} = script_pRFROIAnalysis(conA_data,conB_data,pRFOverlayNames);
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


