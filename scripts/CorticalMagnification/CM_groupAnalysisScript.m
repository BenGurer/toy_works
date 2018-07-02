function CM_groupAnalysisScript

% get study info
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = CM_setupStudyParams;

% use is pc to set data directory - could do in cm_setupStduyparams
% Info.dataDir = '/Volumes/data_PSY/data';
Info.dataDir = 'E:\data';
% E:\data\CorticalMagnification\11108_006
q = char(39);

%% Load subject data
% groupData = struct;
% for iSub = 1:8
for iSub = 5
    % Get subject info
    subjectInfo = get_SubjectInfo_CM(iSub);
    % Subject ID, flatmap names
    saveName = [subjectInfo.subjectID '_data.mat'];
    % move to subject folder
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo.subjectID));
    % load data
    groupData(iSub) = load(saveName);
end

%% Tidy data
% Transform data into tidyVerse
% Row: observation = voxel
% Column: Variable = everything else
% side, group(acquisiton, run), analysis, estimation method, estimate
% for iSub = 1:8
for iSub = 5
    clear data
    data = groupData(iSub).data;
    
    %% Cortical Magnification %%
    pRFanalysisName = ['pRF_', pRFInfo.pRFrestrictROI];
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
    
    subject = repmat(iSub,length(CorticalDistance),1);
    
    T = table(CorticalDistance,Frequency,...
        r2, TuningWidth,...
        Analysis,ROI,...
        subject,...
        'VariableNames',{'CorticalDistance' 'Frequency' 'r2' 'TuningWidth' 'Analysis' 'ROI','Subject'});
    
    writetable(T, [saveName, '_CM.csv'])
    
    %% 7T comparisions %%
    
    %% Concatenated data
    Frequency = [];
    Frequency_kHz = [];
    TuningWidth = [];
    TuningWidth_kHz = [];
    Analysis = [];
    ROI = [];
    r2 = [];
    % voxelPropertyNames: {'Centriod'  'Spread'  'julien_pCF'  'julien_pTW'  'indexMax'}
    % pRFOverlayNames: {'r2'  'PrefCentreFreq'  'rfHalfWidth'}
    estimateFreqNames = {'Centriod', 'julien_pCF', 'indexMax'};
    estimateTuningName = {'Spread', 'julien_pTW', 'NA'};
    
    analysisType = {'GLM', 'pRF'};
    hrfType = {'boxcar', 'Gamma'};
    
    
    pRFanalysisName = ['pRF_', pRFInfo.pRFrestrictROI];
    analysisNames = {'glm_hrfDoubleGamma', 'glm_hrfBoxCar', pRFanalysisName};
    nCons = [32, 8];
    nScans = 4;
    
    for iSide = 1:length(Info.Sides)
        for iGroup = 1:length(glmInfo.groupNames)
            groupName = glmInfo.groupNames{iGroup};
            for iAnal = 1:length(analysisType)
                
                analysisName = analysisNames{iAnal};
                if analysisName ~= analysisNames{3}
                    
                    estimateFreqNames = {'Centriod', 'julien_pCF', 'indexMax'};
                    estimateTuningName = {'Spread', 'julien_pTW', 'NA'};
                    
                else
                    
                    %%%%%%% FIND pRF save names %%%%%%%%%%%%%%%%
                    estimateFreqNames = {'pCF'};
                    estimateTuningName = {'pTW'};
                    
                end                
                
                for iHRF = 1:length(hrfType)                    
                    %                 roiSaveName = [Info.Sides{iSide}, 'GR_' analName{iAnal}];                    
                    roiSaveName = [Info.Sides{iSide}, 'GR_GLM']; % compare using ROI data from derived GLM analysis
                    roiName = [Info.Sides{iSide}, 'GR'];
                    
                    eval(['tempFrequency = data.' roiSaveName '.' groupName '.' analysisName '.' estimateFreqNames{iEst}]);
                    tempFrequency_kHz = funInvNErb(tempFrequency);
                    
                    if ~strcmp('NA',estimateTuningName{iEst})
                        eval(['tempTuningWidth = data.' roiSaveName '.' groupName '.' analysisName '.' estimateTuningName{iEst}]);
                    else
                        tempTuningWidth = nan(length(tempFrequency),1);
                    end
                    tempTuningWidth_kHz = funInvNErb(tempTuningWidth);
                    
                    
                    eval(['tempR2 = data.' roiSaveName '.' groupName '.' analysisName '.r2;']);
                    
                    
                    % seperate loop for scan data
                    if analysisName == analysisNames{1} || analysisNames{2}
                        
                        for iBeta = 1:nCons(1)
                            eval(['beta' mat2str(nCons(1)) ' = data.' roiSaveName '.' groupName '.' analysisName '.betas{iBeta,:}']);
                        end
                        
                    else
                        
                        for iBeta = 1:nCons(1)
                            eval(['betas' mat2str(nCons(1)) ' = nans(length(tempFrequency),1)']);
                        end
                    end
                    
                    nVoxels = length(tempFrequency);
                    
                    tempEstimateFreqName = repmat(estimateFreqNames{iEst},nVoxels,1);
                    tempEstimateTuningName = repmat(estimateTuningName{iEst},nVoxels,1);
                    tempAnalysis = repmat(analysisSaveName{iAnal},nVoxels,1);
                    tempROI = repmat(roiName,nVoxels,1);
                    
                    r2 = [r2; tempR2'];
                    
                    Frequency = [Frequency; tempFrequency'];
                    TuningWidth = [TuningWidth; tempTuningWidth'];
                    Frequency_kHz = [Frequency; tempFrequency_kHz'];
                    TuningWidth_kHz = [TuningWidth; tempTuningWidth_kHz'];
                    
                    if isempty(Analysis)
                        Analysis = tempAnalysis;
                        ROI = tempROI;
                        estimateFreqName = tempEstimateFreqName;
                        estimateTuningName = tempEstimateTuningName;
                    else
                        Analysis = char(Analysis,tempAnalysis);
                        ROI = char(ROI,tempROI);                        
                        estimateFreqName = char(estimateFreqName,tempEstimateFreqName);
                        estimateTuningName = char(estimateTuningName,tempEstimateTuningName);
                    end
                    
                    
                end
                
            end
        end
    end
    
    
    subject = repmat(iSub,length(Frequency),1);
    
    T = table(Frequency,Frequency_kHz,...
        TuningWidth,TuningWidth_kHz,...
        r2, Betas,...
        Analysis,ROI,...
        subject,...
        'VariableNames',{'Frequency' 'Frequency_kHz' 'TuningWidth' 'TuningWidth_kHz' 'r2' 'Betas' 'Analysis' 'ROI','Subject'});
    
    writetable(T, [saveName, '_Comparisions.csv'])
    
    
    %% Scan data
    
    
    
    %% get HRF estimation
    % Two data frames
    % one to average curves
    % one to average fits?
    % OR
    % fit to average curves
    % average and fit in matlab??
    % or export fitted curves as functions and then fit to them?
    
%     data.hrf.x_doubleGamma, data.hrf.x_Gamma, data.hrf.x_dGamma, data.hrf.estimate, data.hrf.deconv, data.hrf.deconvTW
hrfNames = {'Gamma', 'Double Gamma', 'Diff of Gamma', 'BoxCar', 'Deconvolution'};
hrfSaveNames = {'x_Gamma', 'x_doubleGamma', 'x_dGamma', 'deconv'};
nTimepoints = length(data.hrf.estimate.time);

TimePoints = [];
hrfEstimate = [];
hrfEstimateName = [];

get_HRFDoubleGamma(x_doubleGamma,t)
get_HRFGamma(x_Gamma,t)
get_HRFDiffOfGamma(x_dGamma,t)

hrfFunctions = {@get_HRFGamma;
    @get_HRFDoubleGamma;
    @get_HRFDiffOfGamma;
    @get_HRFBoxCar};

for iHRF = 1:length(hrfNames)   
    
    tempTimePoints = data.hrf.estimate.time;
    
    if hrfSaveNames{iHRF} ~= hrfSaveNames{end}
        eval(['temphrfFit = data.hrf.' hrfSaveNames{iHRF} ';']);
        temphrfEstimate = hrfFunctions{iHRF}(temphrfFit,tempTimePoints);
    elseif hrfSaveNames{iHRF} == 'BoxCar'        
        temphrfEstimate = hrfFunctions{iHRF}([2.5, 2.5],tempTimePoints);
    
    else
        eval(['temphrfEstimate = data.hrf.' hrf.SaveNames{iHRF} ';']);        
    end
    
    temphrfEstimateName = repmat(hrfNames{iHRF},nTimepoints,1);
    
    TimePoints = [TimePoints, tempTimePoints];
    hrfEstimate = [hrfEstimate, temphrfEstimate];
    
    if isempty(hrfEstimateName)
        hrfEstimateName = temphrfEstimateName;
    else
        hrfEstimateName = char(hrfEstimateName,temphrfEstimateName);       
    end
       
end
    subject = repmat(iSub,length(TimePoints),1);
    
    T = table(hrfEstimate,TimePoints,...
        hrfEstimateName,...
        subject,...
        'VariableNames',{'hrfEstimate', 'Gamma', 'doubleGamma', 'differenceOfGamma', 'Time(s)', 'Subject'});
    
    writetable(T, [saveName, '_HRF.csv'])
    
    
end

%% get GLM data
% group

% scans


%% get pRF data

%% Estimating HRF
% average HRF estiamte


%% Compare Analysis



%% Compare Acquisiton
% Sparses Vs Continuous
% use best analysis
% Calculate
% Visualize


%% get hrf


end
