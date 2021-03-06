function [ x_doubleGamma, x_Gamma, x_dGamma, estimate, threshold, nVoxels] = script_hrfROIAnalysis(thisView,roiName,glmInfo)
% set to correct group and analysis;
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
roi.scanCoords = getROICoordinates(thisView,roi);
% get ROI estimates
r2data = analysisData.overlays(1).data;
volumeIndices = sub2ind(size(r2data{:}),roi.scanCoords(1,:),roi.scanCoords(2,:),roi.scanCoords(3,:));
[estimate,volumeIndices] = getEstimates(analysisData.d{:}, analysisData.params ,volumeIndices');
% nVoxels = length(volumeIndices);

e = estimate.hdr;
% [index, threshold] = cal_R2threshold(r2data{:}(volumeIndices));
% e = e(:,:,index(:));

% calcculate average HRF
% [ x_doubleGamma, x_Gamma, x_dGamma ] = cal_hrfROIAverage(e,analysisParams)
[ x_doubleGamma, x_Gamma, x_dGamma ] = cal_hrfROIAverage(e,estimate.time,analysisData.d{:});

end
