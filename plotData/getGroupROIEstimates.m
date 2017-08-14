function data = getGroupROIEstimates(thisView,ROI,groupName,analysisName,iScan)
% get data from view
if iScan ~= 0
    thisView = viewSet(thisView,'curGroup',groupName,['curScan=' mat2str(iScan)]);
    else
    thisView = viewSet(thisView,'curGroup',groupName);
    iScan = 1;
end
    thisView = viewSet(thisView,'curAnalysis',viewGet(thisView,'analysisNum',analysisName));
    analysisData{iScan} = viewGet(thisView,'analysis',viewGet(thisView,'analysisNum',analysisName));
    glmData{iScan} = analysisData{iScan}.d{iScan};
    analysisParams{iScan} = analysisData{iScan}.params;
    r2data = analysisData{iScan}.overlays(1).data{iScan};
    % get roi scan coords
    ROI.scanCoords = getROICoordinates(thisView,ROI);
    %get ROI estimates
    volumeIndices = sub2ind(size(r2data),ROI.scanCoords(1,:),ROI.scanCoords(2,:),ROI.scanCoords(3,:));
    [estimate{iScan},volumeIndices] = getEstimates(glmData{iScan} ,analysisParams{iScan} ,volumeIndices');
    nVoxels = length(volumeIndices);
    % save data    
    data.glmData = glmData{iScan};
    data.estimates = estimate{iScan};
    data.betas = squeeze(estimate{iScan}.betas);
    data.betaSte = squeeze(estimate{iScan}.betaSte);
    data.r2 = permute(r2data(volumeIndices),[2 1]);
    
end