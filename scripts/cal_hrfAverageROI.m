function data = cal_hrfAverageROI(e,analysisParams)
    %% Centre max response to 0 and average
    
    hdrTimeMax = 3; % [v i] = max(mean(estimate.hdr(:,:,1),2))
    curve = e.hdr(hdrTimeMax,:,:);
    nOverlays = analysisParams.nhdr;
    resolution = 2;
    overlaySize = size(e.hdr);
    overlaySize = nVoxels;
    %find best frequency using Humphries method
    % remove  values that are less than the average activation across conditions or less than 0
    averageActivation = repmat(mean(curve,3),[1 1 nVoxels]);
    curveNaN=curve;
    curveNaN(curve<averageActivation | curve<0)=NaN;
    %compute the weighted average
    indices = permute(repmat(1:nOverlays,[nVoxels 1]),[3 2 1]);
    bestFrequency = round(resolution*nansum( curveNaN .* indices, 2) ./ nansum(curveNaN,2))/resolution;
    
    % recentre tuning curves using estimated best frequency
    bestFrequency = reshape(bestFrequency,prod(overlaySize),1);
    % curve = reshape(curve,prod(overlaySize),nOverlays);
    uniqueTuningCurveIndices = 1-nOverlays:1/resolution:nOverlays-1;
    tuningCurvesIndices=repmat(permute(uniqueTuningCurveIndices,[1 3 4 2]),[overlaySize 1]);
    tuningCurvesIndices=repmat(permute(uniqueTuningCurveIndices,[1 2]),[overlaySize 1]);
    nTuningCurvesIndices=length(uniqueTuningCurveIndices);
    tuningCurves = nan(1, nTuningCurvesIndices,prod(overlaySize));
    HDRS = nan(analysisParams.nHrfComponents,nTuningCurvesIndices,prod(overlaySize));
    uniqueBestFrequencies = unique(bestFrequency);
    uniqueBestFrequencies=uniqueBestFrequencies(~isnan(uniqueBestFrequencies))';
    for i=uniqueBestFrequencies
        bestFrequencyIndices=find(bestFrequency==i);
        %   tuningCurves(bestFrequencyIndices,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))])=curve(bestFrequencyIndices,:);
        tuningCurves(1,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))],bestFrequencyIndices)=curve(:,:,bestFrequencyIndices);
        HDRS(:,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))],bestFrequencyIndices)=e.hdr(:,:,bestFrequencyIndices);
    end
    % tuningCurves = reshape(tuningCurves,[overlaySize nTuningCurvesIndices]);
    unsmoothedTuningCurves = tuningCurves;
    averageTunindCurves = nansum(unsmoothedTuningCurves,3);
    figure; plot([1-nOverlays:1/resolution:nOverlays-1],averageTunindCurves);
    averageHDRTunindCurves = nansum(HDRS,3);
    figure;
    for i = 1:analysisParams.nHrfComponents
        % plot([1:1:glmData.nHrfComponents],[1-nOverlays:1/resolution:nOverlays-1],averageHDRTunindCurves(i,:));
        plot(averageHDRTunindCurves(i,:));
        legendInfo{i} = ['HDR Component ' num2str(i)];
        hold on
    end
    legend(legendInfo)
    
    figure; waterfall([1-nOverlays:1/resolution:nOverlays-1],[1:1:analysisParams.nHrfComponents],averageHDRTunindCurves);
    figure; surf([1-nOverlays:1/resolution:nOverlays-1],[1:1:analysisParams.nHrfComponents],averageHDRTunindCurves);
    a = permute(averageHDRTunindCurves, [2 1 3]);
    figure; surf([1:1:analysisParams.nHrfComponents],[1-nOverlays:1/resolution:nOverlays-1],a);

end