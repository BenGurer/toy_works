function data = analysisScript_tonotopic(view,concatenationGroupNames,nScans)
% perform tonotopic analysis on fMRI data

% compare:
%   Conditions
%   Between runs
%   Anaylsis type
%       GLM
%       pRF
%       hrfModel

%% Concatenated Group  Analysis
%% GLM analysis
hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
for i = 1:length(concatenationGroupNames)
    for iStim = 1:length(hrfModel)
        analysisName = ['glm_' hrfModel{iStim}];
        thisView = viewSet(thisView,'curGroup',concatenationGroupNames{i});
        [thisView, glmParams] = glmAnalysis(thisView,[],'justGetParams=1','defaultParams=1');
        glmParams.hrfModel = hrfModel{iStim};
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1');
        glmParams.saveName = analysisName;
        glmParams.hrfParams.description = hrfModel{iStim};
        switch hrfModel{iStim}
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
        analysisName = ['glm_' hrfModel{splitHRFmodel} ' - nStim=' mat2str(nStim(iStim)) ' - Scan-' mat2str(iScan)];
        
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
            glmParams.numberContrasts = nStim(iStim);
            glmParams.numberEVs = nStim(iStim);
            glmParams.EVnames = {'1','2','3','4','5','6','7','8'};
        end
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,'justGetParams=1','defaultParams=1',['scanList=' mat2str(iScan)]);
        glmParams.scanParams{iScan}.stimDurationMode = 'From File';
        glmParams.scanParams{iScan}.supersamplingMode =  'Set value';
        glmParams.scanParams{iScan}.designSupersampling = 3;
        glmParams.scanParams{iScan}.acquisitionDelay = .75;
        glmParams.numberContrasts = 0;
        glmParams.parametricTests = 0;
        glmParams.outputEstimatesAsOverlays = 1;
        glmParams.saveName = analysisName;
        [thisView, glmParams] = glmAnalysis(thisView,glmParams,['scanList=' mat2str(iScan)]);
        
        %Tonotopy analysis
        % Index max
        [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str([2:nStim+1])],['scanList=' mat2str(iScan)]);
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
    saveAnalysis(thisView,['glm_' hrfModel{splitHRFmodel} ' - nStim=' mat2str(nStim(iStim)) ' - Scan-' mat2str(iScan)]);
    end
end

% Put back in main function? or above in sub function
% get data from view
data = getGroupROIEstimates(ROI,GroupName,nScan);

% calculate split mean
conditionScans = {[1,3];[2,4]};
for i = 1:length(concatenationGroupNames)
data = cal_splitRun(data(conditionScans{i}(1)),data(conditionScans{i}(2)));
end
% Send to functions to perform analysis
data = plot_ROI_AvB_SplitRuns(run1,run2,run3,run4,Aruns,Bruns,RestrictIndex);
    % contain function to calcuate split mean or make a split mean function and send data onward - prob better

end