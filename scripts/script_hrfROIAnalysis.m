function thisView = script_hrfROIAnalysis(thisView,roiName)
% set to correct group and analysis
% get data from analysis
% use ROI to restrict
% perform ROI analysis - average hrf estimate
% output result

% set view to group we want data from
thisView = viewSet(thisView,'curGroup',glmInfo.groupNames{2});

% get data from analysis
analysisName = 'GLM_Deconv_8bins';
analysisData = get_analysisData(thisView,analysisName);

% get ROI
roi = viewGet(thisView,'roi',roiName);

% restrict data by ROI
% ROIdata = get_ROIdata(analysisData,roi);
    roi.scanCoords = getROICoordinates(thisView,roi);
    %get ROI estimates
    r2data = analysisData.overlays(1).data;
        volumeIndices = sub2ind(size(r2data{:}),roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:));
    [estimate,volumeIndices] = getEstimates(analysisData.d{:}, analysisData.params ,volumeIndices');
    nVoxels = length(volumeIndices);
    
% calcculate average HRF
    data = cal_hrfAverageROI(estimate,analysisData.d{:});
% 
% thisView = viewSet(thisView,'curGroup','MotionComp',['curScan=' mat2str(iScan)]);
%     thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',[analysisSaveName{1} mat2str(iScan)]));
%     analysisData{iScan} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',[analysisSaveName{1} mat2str(iScan)]));
%     glmData{iScan} = analysisData{iScan}.d{iScan};
%     analysisParams{iScan} = analysisData{iScan}.params;
%     r2data = analysisData{iScan}.overlays(1).data{iScan};
%     % get roi scan coords
%     roi{roiNum}.scanCoords = getROICoordinates(thisView,roi{roiNum});
%     %get ROI estimates
%     volumeIndices = sub2ind(size(r2data),roi{roiNum}.scanCoords(1,:),roi{roiNum}.scanCoords(2,:),roi{roiNum}.scanCoords(3,:));
%     [estimate{iScan},volumeIndices] = getEstimates(glmData{iScan} ,analysisParams{iScan} ,volumeIndices');
%     nVoxels = length(volumeIndices);
%     % save data    
%     splitData(iScan).glmData = glmData{iScan};
%     splitData(iScan).estimates = estimate{iScan};
%     splitData(iScan).betas = squeeze(estimate{iScan}.betas);
%     splitData(iScan).betaSte = squeeze(estimate{iScan}.betaSte);
%     splitData(iScan).r2 = permute(r2data(volumeIndices),[2 1]);

end
