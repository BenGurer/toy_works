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
        glmParams.saveName = analysisName;
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
    saveAnalysis(thisView,['glm_' hrfModel{splitHRFmodel} ' - nStim=' mat2str(nStim(iStim)) ' - Scan-' mat2str(iScan)]);
    end
end

% Put back in main function? or above in sub function
% get data from view
% get ROI
ROIName = 'PAC';
pacROI = viewGet(thisView,'roi',ROIName);
% get data from Split runs
% unbinned analysis
for iScan = 1:nScans
ROIEstimatesData_SplitRuns_32Cons{iScan} = getGroupROIEstimates(thisView,pacROI,'MotionComp',['glm_' hrfModel{splitHRFmodel} ' - nStim=32 - Scan-' mat2str(iScan)],iScan);
end
% binned analysis
for iScan = 1:nScans
ROIEstimatesData_SplitRuns_8bins{iScan} = getGroupROIEstimates(thisView,pacROI,'MotionComp',['glm_' hrfModel{splitHRFmodel} ' - nStim=8 - Scan-' mat2str(iScan)],iScan);
end
% get datga from concatenated groups
for iGroup = 1:length(concatenationGroupNames)
ROIEstimatesData_Concat{iGroup} = getGroupROIEstimates(thisView,pacROI,concatenationGroupNames{iGroup},['glm_' hrfModel{2}],0);
end

% calculate split mean
conditionScans = {[1,3];[2,4]};

conditionScans = {[2,4];[1,3]};
for iGroup = 1:length(concatenationGroupNames)
runA_32cons{iGroup} = ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(1)}.betas;
runB_32cons{iGroup} = ROIEstimatesData_SplitRuns_32Cons{conditionScans{iGroup}(2)}.betas;

runA_8bins{iGroup} = ROIEstimatesData_SplitRuns_8bins{conditionScans{iGroup}(1)}.betas;
runB_8bins{iGroup} = ROIEstimatesData_SplitRuns_8bins{conditionScans{iGroup}(2)}.betas;

restrictIndex = true(1,length(runA_32cons{iGroup}));

% replace voxels outside of restriction with nans
runA_32cons{iGroup}(~restrictIndex) = nan;
runB_32cons{iGroup}(~restrictIndex) = nan;
runA_8bins{iGroup}(~restrictIndex) = nan;
runB_8bins{iGroup}(~restrictIndex) = nan;

runA_mv{iGroup} = cal_movingAverage(runA_32cons{iGroup});
runB_mv{iGroup} = cal_movingAverage(runB_32cons{iGroup});

runA_Peak_mv{iGroup} = cal_voxel_properties(runA_mv{iGroup});
runB_Peak_mv{iGroup} = cal_voxel_properties(runB_mv{iGroup});

roi_av{iGroup} = cal_ROI_pTW_av(runA_32cons{iGroup},runB_32cons{iGroup});
roi_av_mv{iGroup} = cal_ROI_pTW_av(runA_mv{iGroup},runB_mv{iGroup});
roi_av_bin{iGroup} = cal_ROI_pTW_av(runA_8bins{iGroup},runB_8bins{iGroup});

% [condition_splitMean, condition_splitMean_peak, condition_Cntrd, condition_Sprd, condition_pRF, voxel_Mean_Peak, Voxel_Cntrd, Voxel_Sprd, Voxel_pRF] = cal_splitMean(betas_A,betas_B) 
% [condition_splitMean_mv{iGroup}, condition_splitMean_peak_mv{iGroup}, condition_Cntrd_mv{iGroup}, condition_Sprd_mv{iGroup}, condition_pRF_mv{iGroup}, voxel_Mean_Peak_mv{iGroup}, Voxel_Cntrd_mv{iGroup}, Voxel_Sprd_mv{iGroup}, Voxel_pRF_mv{iGroup}] = cal_splitMean(runA_mv{iGroup},runB_mv{iGroup});
% [condition_splitMean_bin{iGroup}, condition_splitMean_peak_bin{iGroup}, condition_Cntrd_bin{iGroup}, condition_Sprd_bin{iGroup}, condition_pRF_bin{iGroup}, voxel_Mean_Peak_bin{iGroup}, Voxel_Cntrd_bin{iGroup}, Voxel_Sprd_bin{iGroup}, Voxel_pRF_bin{iGroup}] = cal_splitMean(runA_8bins{iGroup},runB_8bins{iGroup});
[condition_splitMean_mv{iGroup}, ROI_data_mv{iGroup}, Voxel_data_mv{iGroup}] = cal_splitMean(runA_mv{iGroup},runB_mv{iGroup});
 [condition_splitMean_bin{iGroup}, ROI_data_bin{iGroup}, Voxel_data_bin{iGroup}] = cal_splitMean(runA_8bins{iGroup},runB_8bins{iGroup});

end

% %% Plot results using:
% plot_compareConditions_pCF
% plot_compareConditions_pTW
% plot_compareConditions_ROIav

%% Compare average ROI beta weights
% ADD ste to av ROI plots
% recentre AvROI plot
plot_compareConditions_ROIav(roi_av_bin{1},roi_av_bin{2})
plot_compareConditions_ROIav(roi_av_mv{1},roi_av_mv{2})
plot_compareConditions_ROIav(roi_av{1},roi_av{2})

% % Condition A TRUE compare conditions
% ConB_runA = ROIEstimatesData_SplitRuns_8bins{conditionScans{1}(1)}.betas;
% ConB_runB = ROIEstimatesData_SplitRuns_8bins{conditionScans{1}(2)}.betas;

% plot_compareCondtions_byROIAverage_pTW_ConditionATRUE(condition_splitMean_bin{1},condition_splitMean_peak_bin{1},runA_8bins,runB_8bins)
% [conA, conB] = cal_compareConditions_byVoxel_pTW(voxel_Mean_bin{1},voxel_Mean_bin{2});

%% Compare ROI pCF
% Compare conditions to runs
% compare runs
plot_compareConditions_pCF(runA_Peak_mv{1}, runB_Peak_mv{1},24);
plot_compareConditions_pCF(runA_Peak_mv{2}, runB_Peak_mv{2},24);

%% need to update below to work on structures%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare conditions
plot_compareConditions_pCF(voxel_Mean_Peak_bin{1},voxel_Mean_Peak_bin{2},8);
plot_compareConditions_pCF(Voxel_Cntrd_bin{1},Voxel_Cntrd_bin{2},8);

plot_compareConditions_pCF(voxel_Mean_Peak_mv{1},voxel_Mean_Peak_mv{2},28);
plot_compareConditions_pCF(Voxel_Cntrd_mv{1},Voxel_Cntrd_mv{2},28);

% plot pRFs and calculate centriod
plot_compareConditions_pCF(condition_pRF_mv{1},condition_pRF_mv{2},28);
plot_compareConditions_pCF(Voxel_pRF_mv{1},Voxel_pRF_mv{2},28);
plot_compareConditions_pCF(condition_pRF_bin{1},condition_pRF_bin{2},8);
plot_compareConditions_pCF(Voxel_pRF_bin{1},Voxel_pRF_bin{2},8);

%% Compare ROI pTW
plot_compareConditions_pTW(condition_splitMean_bin{1},condition_splitMean_bin{2},2);
plot_compareConditions_pTW(condition_splitMean_mv{1},condition_splitMean_mv{2},4);

%% Compare ROI pTW assuming condition A is TRUE
ConBROIpTW_bin = cal_ConBROIpTW_ConAVoxelIndex(voxel_Mean_Peak_bin{1},runA_8bins{2},runB_8bins{2});
plot_compareConditions_pTW(condition_splitMean_bin{1}, ConBROIpTW_bin,2);

ConBROIpTW_mv = cal_ConBROIpTW_ConAVoxelIndex(voxel_Mean_Peak_mv{1},runA_mv{2},runB_mv{2});
plot_compareConditions_pTW(condition_splitMean_mv{1}, ConBROIpTW_mv,4);

% plot_compareConditions_byROIAverage(conA, conB,2);

% Send to functions to perform analysis
% data = plot_ROI_AvB_SplitRuns(run1,run2,run3,run4,Aruns,Bruns,RestrictIndex);
    % contain function to calcuate split mean or make a split mean function and send data onward - prob better
    
    % reshape and plot as overlay? see most effected areas on the cortex

end