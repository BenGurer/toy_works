function data = getROIdata_GLM(thisView,analysisData,ROI,iScan)
    %
    %   usage: getROIdata_GLM(thisView,analysisData,ROI,iScan)
    %      by: Ben Gurer
    %    date: 05/10/2017
    % purpose: Get GLM data within ROI
    %   input: mrTools view, analysis data, ROI, scan number (if defined)
    %
if iScan ~= 0
else
    iScan = 1;
end
% use r2 data to set size of data set for indexing
r2data = analysisData.overlays(1).data{iScan};
% get roi scan coords
ROI.scanCoords = getROICoordinates(thisView,ROI);
% get ROI estimates
volumeIndices = sub2ind(size(r2data),ROI.scanCoords(1,:),ROI.scanCoords(2,:),ROI.scanCoords(3,:));
[estimates,volumeIndices] = getEstimates(analysisData.d{iScan} , analysisData.params, volumeIndices');
% save data
data.glmData = analysisData.d{iScan};
data.estimates = estimates;
data.betas = squeeze(estimates.betas);
data.betaSte = squeeze(estimates.betaSte);
data.r2 = permute(r2data(volumeIndices),[2 1]);
end