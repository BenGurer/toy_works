function data = getROIdata(thisView,analysisData,ROI,iScan)
if iScan ~= 0
else
    iScan = 1;
end
% glmData{iScan} = analysisData{iScan}.d{iScan};
% analysisParams{iScan} = analysisData.params;
r2data = analysisData.overlays(1).data{iScan};
% get roi scan coords
ROI.scanCoords = getROICoordinates(thisView,ROI);
%get ROI estimates
volumeIndices = sub2ind(size(r2data),ROI.scanCoords(1,:),ROI.scanCoords(2,:),ROI.scanCoords(3,:));
[estimates,volumeIndices] = getEstimates(analysisData.d{iScan} , analysisData.params, volumeIndices');
% nVoxels = length(volumeIndices);
% save data
data.glmData = analysisData.d{iScan};
data.estimates = estimates;
data.betas = squeeze(estimates.betas);
data.betaSte = squeeze(estimates.betaSte);
data.r2 = permute(r2data(volumeIndices),[2 1]);
end