function data = analysisScript_tonotopic(thisView,concatenationGroupNames,ROIName,nScans)
% perform tonotopic analysis on fMRI data

% compare:
%   Conditions
%   Between runs
%   Anaylsis type
%       GLM
%       pRF
%       hrfModel

%% Concatenated Group  Analysis
%
% thisView = getMLRView;
% nScans = 4;
%
% ROIName = 'AC';
% concatenationGroupNames = {'ConcatenationSparse', 'ConcatenationCont'};
% concatenationGroupNames = {'ConcatenationHLsim', 'ConcatenationNH'};

%% GLM analysis
hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
for i = 1:length(concatenationGroupNames)
    for iHRF = 1:length(hrfModel)
        analysisName = ['glm_' hrfModel{iHRF}];
        thisView = viewSet(thisView,'curGroup',concatenationGroupNames{i});
        [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
        glmParams.hrfModel = hrfModel{iHRF};
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
        glmParams.saveName = analysisName;
        glmParams.hrfParams.description = hrfModel{iHRF};
        switch hrfModel{iHRF}
            case 'hrfBoxcar'
                glmParams.hrfParams.delayS =  2.5;
                glmParams.hrfParams.durationS = 2.5;
            case 'hrfDoubleGamma'
                glmParams.hrfParams.x =  4;
                glmParams.hrfParams.y = 11;
                glmParams.hrfParams.z =  4;
        end
        glmParams.scanParams{1}.stimDurationMode = 'From file';
        glmParams.scanParams{1}.supersamplingMode =  'Set value';
        glmParams.scanParams{1}.designSupersampling = 3;
        glmParams.scanParams{1}.acquisitionDelay = .75;
        glmParams.computeTtests = 1;
        glmParams.numberContrasts  = 1;
        glmParams.contrasts = ones(1,32);
        glmParams.numberFtests  = 1;
        glmParams.fTestNames{1, 1} = 'fTest - all conditions';
        glmParams.restrictions{1, 1} = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0;...
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
        glmParams.alphaContrastOverlay = 'Uncorrected';
        glmParams.parametricTests = 1;
        glmParams.fweAdjustment = 0;
        glmParams.fdrAdjustment = 0;
        glmParams.outputStatistic = 0;
        glmParams.numberContrasts = 0;
        glmParams.outputEstimatesAsOverlays = 1;
        [thisView, glmParams] = glmAnalysis(thisView,glmParams);
        
        %Tonotopy analysis
        % Index max
        [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:33])]);
        params.combineFunction='indexMax';
        params.nOutputOverlays=2;
        [thisView,params] = combineTransformOverlays(thisView,params);
        curOverlay=viewGet(thisView,'curOverlay');
        thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);
        
        % Weighted mean and corrected weighted mean
        params.combineFunction='weightedMeanStd';
        params.nOutputOverlays=4;
        [thisView,params] = combineTransformOverlays(thisView,params);
        curOverlay=viewGet(thisView,'curOverlay');
        thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
        thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
        thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
        thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);
        
        %% save analysis
        saveAnalysis(thisView,analysisName)
    end
end

%% Split Run Analysis
%% GLM Analysis
nStim = [8 32];
% analysisSaveName = {'GLM Box Car - ALL CONS - Scan ','pRF - ALL CONS - Scan '};
% Set hrf type based on acquisition type
splitHRFmodel = 2;
for iScan = 1:nScans
    for iStim = 1:length(nStim)
        analysisName_split{iScan} = ['glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '_Scan_' mat2str(iScan)];
        
        thisView = viewSet(thisView,'curGroup','MotionComp',['curScan=' mat2str(iScan)]);
        [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
        glmParams.hrfModel = hrfModel{splitHRFmodel};
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
        %         glmParams.saveName = analysisName;
        glmParams.hrfParams.description = hrfModel{splitHRFmodel};
        switch hrfModel{splitHRFmodel}
            case 'hrfBoxcar'
                glmParams.hrfParams.delayS =  2.5;
                glmParams.hrfParams.durationS = 2.5;
            case 'hrfDoubleGamma'
                glmParams.hrfParams.x =  4;
                glmParams.hrfParams.y = 11;
                glmParams.hrfParams.z =  4;
        end
        
        if nStim(iStim) == 8
            glmParams.scanParams{1, iScan}.preprocess  = 'binStimFreq';
            glmParams.EVnames = {'1','2','3','4','5','6','7','8'};
        end
        glmParams.numberContrasts = nStim(iStim);
        glmParams.numberEVs = nStim(iStim);
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
        glmParams.scanParams{iScan}.stimDurationMode = 'From File';
        glmParams.scanParams{iScan}.supersamplingMode =  'Set value';
        glmParams.scanParams{iScan}.designSupersampling = 3;
        glmParams.scanParams{iScan}.acquisitionDelay = .75;
        glmParams.numberContrasts = 0;
        glmParams.parametricTests = 0;
        glmParams.outputEstimatesAsOverlays = 1;
        glmParams.saveName = analysisName_split{iScan};
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,['scanList=' mat2str(iScan)]);
        
        %Tonotopy analysis
        % Index max
        [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:nStim(iStim)+1])],['scanList=' mat2str(iScan)]);
        params.combineFunction='indexMax';
        params.nOutputOverlays=2;
        [thisView,params] = combineTransformOverlays(thisView,params);
        curOverlay=viewGet(thisView,'curOverlay');
        thisView = viewSet(thisView,'overlaycolorrange',[0 nStim(iStim)],curOverlay-1);
        
        % Weighted mean - if all stimuli
        if nStim(iStim) == 32
            params.combineFunction='weightedMeanStd';
            params.nOutputOverlays=4;
            [thisView,params] = combineTransformOverlays(thisView,params);
            curOverlay=viewGet(thisView,'curOverlay');
            thisView = viewSet(thisView,'overlaycolorrange',[0 nStim(iStim)],curOverlay-3);
            thisView = viewSet(thisView,'overlaycolorrange',[0 nStim(iStim)],curOverlay-2);
            thisView = viewSet(thisView,'overlaycolorrange',[0 nStim(iStim)*1.25],curOverlay-1);
            thisView = viewSet(thisView,'overlaycolorrange',[0 nStim(iStim)*1.25],curOverlay);
        end
    end
end
% save analysis
for iScan = 1:nScans
    for iStim = 1:length(nStim)
        saveAnalysis(thisView,analysisName_split{iScan});
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% change to get all data then get from ROI %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put back in main function? or above in sub function
% get data from view
% get ROI
% ROIName = 'PAC';

% Change to get all data then restict by ROI
thisView = getMLRView;
ROIName = 'RightAC';
ROIName = 'LeftAC';
ROIName = 'pAC';
ROIName = 'RIGHT';
ROIName = 'AntHGAC';
pacROI = viewGet(thisView,'roi',ROIName);
% get data from Split runs
% unbinned analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% THIS NEEDS TO BE IMPROVED - GET DATA ONCE! %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load Group
% 2. Get scan data - loop here
% 3. Restrict by ROI coords

% for iScan = 1:nScans
%     ROIEstimatesData_SplitRuns_32Cons{iScan} = getGroupROIEstimates(thisView,pacROI,'MotionComp',['glm_' hrfModel{splitHRFmodel} '_nCons_32_Scan_' mat2str(iScan)],iScan);
% end
% % binned analysis
% for iScan = 1:nScans
%     ROIEstimatesData_SplitRuns_8bins{iScan} = getGroupROIEstimates(thisView,pacROI,'MotionComp',['glm_' hrfModel{splitHRFmodel} '_nCons_8_Scan_' mat2str(iScan)],iScan);
% end
% % get datga from concatenated groups
% for iGroup = 1:length(concatenationGroupNames)
%     ROIEstimatesData_Concat{iGroup} = getGroupROIEstimates(thisView,pacROI,concatenationGroupNames{iGroup},['glm_' hrfModel{2}],0);
% end




%% Save data in structures labelled by roi
%% ROI names
% Create ROIs with these names:
% RightAC
% RightPosAC
% RightAntAc
% LeftAC
% LeftPosAC
% LeftAntAC
% AC
ROInames = {'RightAC','RightPosAC','RightAntAC','LeftAC','LeftPosAC','LeftAntAC','AC'};

ROInames = {'RightAC','RightPosAC'};

q = char(39);
for iROI = 1:length(ROInames)
eval([ROInames{iROI} ' = struct;']);
eval([ROInames{iROI} '.roi = viewGet(thisView,' q 'roi' q ',ROInames{iROI});']);
end
% function to create figure based on paper size - poster or paper
% send figure handle to plotting functions - can be used with subplot and plot
% create strucutre for rois or containing rois

% following analysis - save roi results underone structure of subject ID
% plots I want
% ROI average beta weight
% ROI TW conA v B

%% SHOULD SET ANALYSIS NAME ONCE AND PASS IT ON
% loop - create analysis names
% loop length of analysis names
% needs to name based on group - one structure

%% Get data from individual scans
thisView = viewSet(thisView,'curGroup','MotionComp');
for iStim = 1:length(nStim)
    for iScan = 1:nScans
        analysisName = ['glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '_Scan_' mat2str(iScan)];
        eval(['GLManalysisData_nCons_' mat2str(nStim(iStim)) '{iScan} = getGroupAnalysisData(thisView,analysisName);'])
    end
end
% get data from ROIs
for iROI = 1:length(ROInames)
    for iScan = 1:nScans
        for iStim = 1:length(nStim)
            eval([ROInames{iROI} '.glm_' hrfModel{splitHRFmodel} '_nCons_' mat2str(nStim(iStim)) '{iScan} = getROIdata(thisView,GLManalysisData_nCons_' mat2str(nStim(iStim)) '{iScan},' ROInames{iROI} '.roi,iScan);']);
        end
    end
end

% now get concat group data
thisView = viewSet(thisView,'curGroup',groupName);

% get datga from concatenated groups
for iGroup = 1:length(concatenationGroupNames)
    
thisView = viewSet(thisView,'curGroup',concatenationGroupNames{iGroup});

    ROIEstimatesData_Concat{iGroup} = getGroupROIEstimates(thisView,pacROI,concatenationGroupNames{iGroup},['glm_' hrfModel{2}],0);
end
%%
% get experimental data

%% now pass data to function to plot
data = roiDataAnalysis(ROIdata);
% loop over ROIs

% 
% ROIEstimatesData_SplitRuns_32Cons{iScan} = getROIdata(thisView,GLManalysisData_32cons{iScan},pacROI,iScan);
%     ROIEstimatesData_SplitRuns_8bins{iScan} = getROIdata(thisView,GLManalysisData_8bins{iScan},pacROI,iScan);

% change to plot equally spaced by stim condition but lablled by kHz - ie ERB Scaling

% decide which plots I need
% science and explaining
% explaining - mean vs split mean

% get GLM analysis data - 32 conditions
%     analysisName = ['glm_' hrfModel{splitHRFmodel} '_nCons_32_Scan_' mat2str(iScan)];
%     GLManalysisData_32cons{iScan} = getGroupAnalysisData(thisView,analysisName,iScan);
%     % get GLM analysis data - 8 bins
%     analysisName = ['glm_' hrfModel{splitHRFmodel} '_nCons_8_Scan_' mat2str(iScan)];
%     GLManalysisData_8bins{iScan} = getGroupAnalysisData(thisView,analysisName,iScan);
%

% ROI names
% RightAC
% RightPosAC
% RightAntAc
% LeftAC
% LeftPosAC
% LeftAntAC
% AC

% calculate split mean
% 7 T
% conditionScans = {[1,3];[2,4]};
% 3 T
conditionScans = {[2,4];[1,3]};
conditionScans = {[6,8];[5,7]};

restrictIndex = true(1,length(ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(1)}.betas));
r2Threshold = (max(ROIEstimatesData_Concat{2}.r2)-min(ROIEstimatesData_Concat{2}.r2)).*0.7;
restrictIndex = ROIEstimatesData_Concat{2}.r2 > r2Threshold;
nVoxels = sum(restrictIndex);

% replace voxels outside of restriction with nans
% runA_32cons{iGroup}(~restrictIndex) = nan;
% runB_32cons{iGroup}(~restrictIndex) = nan;
% runA_8bins{iGroup}(~restrictIndex) = nan;
% runB_8bins{iGroup}(~restrictIndex) = nan;

runA_32cons = [];
runB_32cons = [];
runA_8bins = [];
runB_8bins = [];
runA_32cons_restricted = [];
runB_32cons_restricted = [];
runA_8bins_restricted = [];
runB_8bins_restricted = [];

runA_mv = [];
runB_mv = [];
for iGroup = 1:length(concatenationGroupNames)
    runA_32cons{iGroup} = ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(1)}.betas;
    runB_32cons{iGroup} = ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(2)}.betas;
    
    runA_32cons{iGroup}(isnan(runA_32cons{iGroup})) = [];
    runB_32cons{iGroup}(isnan(runB_32cons{iGroup})) = [];
    
    runA_8bins{iGroup} = ROIEstimatesData_SplitRuns_8bins{conditionScans{iGroup}(1)}.betas;
    runB_8bins{iGroup} = ROIEstimatesData_SplitRuns_8bins{conditionScans{iGroup}(2)}.betas;
    
    runA_32cons_restricted{iGroup} = runA_32cons{iGroup}(:,restrictIndex);
    runB_32cons_restricted{iGroup} = runB_32cons{iGroup}(:,restrictIndex);
    runA_8bins_restricted{iGroup} = runA_8bins{iGroup}(:,restrictIndex);
    runB_8bins_restricted{iGroup} = runB_8bins{iGroup}(:,restrictIndex);
    
    runA_mv{iGroup} = cal_movingAverage(runA_32cons_restricted{iGroup});
    runB_mv{iGroup} = cal_movingAverage(runB_32cons_restricted{iGroup});
    
    runA_Peak_mv{iGroup} = cal_voxel_properties(runA_mv{iGroup});
    runB_Peak_mv{iGroup} = cal_voxel_properties(runB_mv{iGroup});
    
    roi_av{iGroup} = cal_ROI_pTW_av(runA_32cons_restricted{iGroup},runB_32cons_restricted{iGroup});
    roi_av_mv{iGroup} = cal_ROI_pTW_av(runA_mv{iGroup},runB_mv{iGroup});
    roi_av_bin{iGroup} = cal_ROI_pTW_av(runA_8bins_restricted{iGroup},runB_8bins_restricted{iGroup});
    
    % [condition_splitMean{iGroup}, ROI_data{iGroup}, Voxel_data{iGroup}, totalROIpCF{iGroup}] = cal_splitMean(runA_32cons_restricted{iGroup},runB_32cons_restricted{iGroup});
    [condition_splitMean_mv{iGroup}, ROI_data_mv{iGroup}, Voxel_data_mv{iGroup}, totalROIpCF_mv{iGroup}] = cal_splitMean(runA_mv{iGroup},runB_mv{iGroup});
    [condition_splitMean_bin{iGroup}, ROI_data_bin{iGroup}, Voxel_data_bin{iGroup}, totalROIpCF_bin{iGroup}] = cal_splitMean(runA_8bins_restricted{iGroup},runB_8bins_restricted{iGroup});
end

%% get stimulus properties
lowFreqkHz = 0.1;
highFreqkHz = 8;
nStim = 32;
[stimFreqs, stimFreqs_bin, stimFreqs_mv] = convertStimIDtoFrequency(lowFreqkHz,highFreqkHz,nStim);

[stimLevel_SL, maskingLevel] = calStimulusSensationLevel(stimFreqs);

figure('color',[1 1 1]); plot(stimFreqs,stimLevel_SL);xlabel('Frequency (kHz)'); ylabel('Sensation Level (dB SL)')

stimLevel_SL_bin = [];
binSize = 4;
c = 1;
for i = 1:length(stimLevel_SL)/binSize
    stimLevel_SL_bin(i) = mean(stimLevel_SL(c:c + (binSize-1)));
    c = c + binSize;
end

% moving average

nBins = 8;
windowAvSize = length(stimLevel_SL)/nBins;

if isreal(windowAvSize) && rem(windowAvSize,1)==0
    
    loopLength = length(stimLevel_SL) - windowAvSize;
    stimLevel_SL_mv = zeros(1,loopLength);
    
    for i = 1:loopLength
        stimLevel_SL_mv(i) = mean(stimLevel_SL(i:i+windowAvSize-1));
    end
    
else
    error('Moving average window not an integer')
end

% %% Plot results using:
% plot_compareConditions_pCF
% plot_compareConditions_pTW
% plot_compareConditions_ROIav

%% Compare average ROI beta weights
% ADD ste to av ROI plots
% recentre AvROI plot
plot_compareConditions_ROIav(roi_av_bin{1},roi_av_bin{2},stimFreqs_bin)
plot_compareConditions_ROIav(roi_av_mv{1},roi_av_mv{2},stimFreqs_mv)
plot_compareConditions_ROIav(roi_av{1},roi_av{2},stimFreqs)

%% Compare ROI pCF
% Compare conditions to runs
% compare runs
plot_compareConditions_pCF(runA_Peak_mv{1}, runB_Peak_mv{1},24);
plot_compareConditions_pCF(runA_Peak_mv{2}, runB_Peak_mv{2},24);

% compare conditions
plot_compareConditions_pCF(Voxel_data_bin{1}.Cntrd,Voxel_data_bin{2}.Cntrd,8);
plot_compareConditions_pCF(Voxel_data_mv{1}.Cntrd,Voxel_data_mv{2}.Cntrd,28);

% plot pRFs and calculate centriod
plot_compareConditions_pCF(Voxel_data_bin{1}.pRF_params(:,1)',Voxel_data_bin{2}.pRF_params(:,1)',8);
plot_compareConditions_pCF(Voxel_data_mv{1}.pRF_params(:,1)',Voxel_data_mv{2}.pRF_params(:,1)',28);

plot_compareConditions_pCF(Voxel_data_bin{1}.pRF_cntrd,Voxel_data_bin{2}.pRF_cntrd,8);
plot_compareConditions_pCF(Voxel_data_mv{1}.pRF_cntrd,Voxel_data_mv{2}.pRF_cntrd,28);

%% Compare ROI pTW
plot_compareConditions_pTW(condition_splitMean_bin{1},condition_splitMean_bin{2},2);
plot_compareConditions_pTW(condition_splitMean_mv{1},condition_splitMean_mv{2},4);

%% Compare ROI pTW assuming condition A is TRUE
ConBROIpTW_bin = cal_ConBROIpTW_ConAVoxelIndex(Voxel_data_bin{1}.Mean_Peak,runA_8bins{2},runB_8bins{2});
plot_compareConditions_pTW(condition_splitMean_bin{1}, ConBROIpTW_bin,2,stimFreqs_bin);

ConBROIpTW_mv = cal_ConBROIpTW_ConAVoxelIndex(Voxel_data_mv{1}.Mean_Peak,runA_mv{2},runB_mv{2});
plot_compareConditions_pTW(condition_splitMean_mv{1}, ConBROIpTW_mv,4,stimFreqs_mv);


%% Compare ROI pTW peak weights



% [peakWeightDiff_bin, peakWeightRatio_bin] = cal_pTWpeakWeight(condition_splitMean_bin{1}, ConBROIpTW_bin);
% figure
% plot(peakWeightDiff_bin);
% [peakWeightDiff_mv, peakWeightRatio_mv] = cal_pTWpeakWeight(condition_splitMean_mv{1}, ConBROIpTW_mv);
% figure
% plot(peakWeightDiff_mv);



peakWeightRatio = roi_av{2}./roi_av{1};
peakWeightRatio_bin = roi_av_bin{2}./roi_av_bin{1};
peakWeightRatio_mv = roi_av_mv{2}./roi_av_mv{1};

ratio2Plot = [peakWeightRatio(stimLevel_SL~=50); mean(peakWeightRatio(stimLevel_SL==50))];

level2Plot = [stimLevel_SL(stimLevel_SL~=50) 50];


ratio2Plot_bin = [peakWeightRatio_bin(stimLevel_SL_bin~=50); mean(peakWeightRatio_bin(stimLevel_SL_bin==50))];
level2Plot_bin = [stimLevel_SL_bin(stimLevel_SL_bin~=50) 50];


ratio2Plot_mv  = [peakWeightRatio_mv(stimLevel_SL_mv~=50); mean(peakWeightRatio_mv(stimLevel_SL_mv==50))];
level2Plot_mv  = [stimLevel_SL_mv(stimLevel_SL_mv ~=50) 50];
% stimFreqs2Plot_mv = [stimFreqs_mv(stimLevel_SL_mv~=50); mean(stimFreqs_mv(stimLevel_SL_mv==50))];

figure; scatter(level2Plot,ratio2Plot)
figure; scatter(level2Plot_bin,ratio2Plot_bin)


%% dB SL vs beta weight
figure('color',[1 1 1]); scatter(level2Plot_mv ,ratio2Plot_mv )
hold on
fit = polyfit(level2Plot_mv',ratio2Plot_mv,1);
plot(level2Plot_mv,polyval(fit,level2Plot_mv));
% correlation = corrcoef([level2Plot_mv' ratio2Plot_mv]);
% text(25,0.8,sprintf('Correlation = %.2f \n y = %.2fx + %.2f',correlation(2), fit(1), fit(2)))
% plot(f,f,'k--')
errorbar(mean(stimLevel_SL(stimLevel_SL==50)),mean(peakWeightRatio(stimLevel_SL==50)),std(peakWeightRatio(stimLevel_SL==50)))
xlabel('Stimulus Sensation Level (dB SL)'); ylabel('Average Beta Weight Ratio (ConB / ConA)');



figure; scatter(stimLevel_SL,peakWeightRatio)
figure; scatter(stimLevel_SL_bin,peakWeightRatio_bin)
figure; scatter(stimLevel_SL_mv,peakWeightRatio_mv)

figure; scatter(stimLevel_SL_bin,peakWeightDiff_bin)
figure; scatter(stimLevel_SL_mv,peakWeightDiff_mv)
%
% MaskingLevelatStimuli = getMaskingLevelatFrequency();
% plot_diffPeakWeightVsNoiseLevel(peakWeightDiff,MaskingLevelatStimuli);

%% TO DO
% non split mean to justify split mean

% % create function to convert to frequency
% convert stimuli to frequency - invNERB plus nERB(0.1kHZ)
% convert to binned - 8 groups - 4 per - middle between bin lims
% figure out moving average
% use this for plotting
% get masking level at frequecny

% reshape and plot as overlay? see most effected areas on the cortex
% threshold based on statisitcal test

end