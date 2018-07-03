function [ x_doubleGamma, x_Gamma, x_dGamma, hrf_Deconv, HRF_TW_est] = cal_hrfROIAverage(e,t,analysisParams)
    %
    %   usage: cal_hrfROIAverage
    %      by: Ben Gurer
    %    date: 03/28/2018
    % purpose: calculate average HRF within given ROI
    %   input: voxel estimates from ROI, analysis parameters
    %  output: Fitted HRF params for double gamma, gamma and difference of
    %  gamma models. ROI Average HRF, ROI average Tuning width and HRF
    %  surface
    %
%% 
%% Centre max response to 0 and average
%     max(a,[],2)
% take mean of all and find max time point and use that
%     hdrMaxTimePoint = 3; % [v i] = max(mean(estimate.hdr(:,:,1),2))
[v hdrMaxTimePoint] = max(max(mean(e,3)'));
curve = e(hdrMaxTimePoint,:,:);
nOverlays = analysisParams.nhdr;
resolution = 1; %resolution = 2;
overlaySize = size(e);
% overlaySize = nVoxels;
nVoxels = overlaySize(3);
%find best frequency using Humphries method
% remove  values that are less than the average activation across conditions or less than 0
averageActivation = repmat(mean(curve,3),[1 1 nVoxels]);
curveNaN=curve;
curveNaN(curve<averageActivation | curve<0)=NaN;
%compute the weighted average
indices = permute(repmat(1:nOverlays,[nVoxels 1]),[3 2 1]);
bestFrequency = round(resolution*nansum( curveNaN .* indices, 2) ./ nansum(curveNaN,2))/resolution;

% recentre tuning curves using estimated best frequency
bestFrequency = reshape(bestFrequency,prod(nVoxels),1);
% curve = reshape(curve,prod(overlaySize),nOverlays);
uniqueTuningCurveIndices = 1-nOverlays:1/resolution:nOverlays-1;
tuningCurvesIndices=repmat(permute(uniqueTuningCurveIndices,[1 3 4 2]),[nVoxels 1]);
tuningCurvesIndices=repmat(permute(uniqueTuningCurveIndices,[1 2]),[nVoxels 1]);
nTuningCurvesIndices=length(uniqueTuningCurveIndices);
tuningCurves = nan(1, nTuningCurvesIndices,prod(nVoxels));
HDRS = nan(analysisParams.nHrfComponents,nTuningCurvesIndices,prod(nVoxels));
uniqueBestFrequencies = unique(bestFrequency);
uniqueBestFrequencies=uniqueBestFrequencies(~isnan(uniqueBestFrequencies))';
for i=uniqueBestFrequencies
    bestFrequencyIndices=find(bestFrequency==i);
    %   tuningCurves(bestFrequencyIndices,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))])=curve(bestFrequencyIndices,:);
    tuningCurves(1,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))],bestFrequencyIndices)=curve(:,:,bestFrequencyIndices);
    HDRS(:,find(abs(uniqueTuningCurveIndices-1+i)<1e-6)+[0 cumsum(repmat(resolution,1,nOverlays-1))],bestFrequencyIndices)=e(:,:,bestFrequencyIndices);
end
% tuningCurves = reshape(tuningCurves,[overlaySize nTuningCurvesIndices]);
% unsmoothedTuningCurves = tuningCurves;
averageTunindCurves = nansum(tuningCurves,3);
averageHDRTunindCurves = nansum(HDRS,3);

%% Plotting fitting of HRF or HRF deconvolution estimate
% get data ready/ in the right format
timePoints = t;
freqBin = 1-nOverlays : 1/resolution : nOverlays-1;
% permute "average hrf tunding curve" matrix dimentions (from frequency,time to time, frequency)
HRF_TW_est = permute(averageHDRTunindCurves, [2 1 3]);
% normalise to allow analysis to scale and offset
hrf_Deconv = HRF_TW_est(7,:)/max(HRF_TW_est(7,:));


% Now plot
figure
subplot(3,2,1)
waterfall(t,freqBin,HRF_TW_est)
% % surf([1:1:analysisParams.nHrfComponents],[1-nOverlays:1/resolution:nOverlays-1],HRF_TW_est);
% 
% p_fmribHRF = [6 12 0.9 0.9 0.35 1]; %guess
% opts = optimset('MaxFunEvals', 500, 'Display', 'off');
% [x_fmribHRF, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFfmrib, p_fmribHRF, t, hrf_Deconv, [], [], opts);
% % disp(x_fmribHRF);
% subplot(3,2,2)
% fittedHRF = get_HRFfmrib(x_fmribHRF,t);
% plot(t,fittedHRF)
% hold on
% plot(t,hrf_Deconv)
% title('fmribHRF')
% legend('fitted','deconv')


%% GLM Double Gamma
% gammaAshapeA = x(1); % must be positive
% gammaBshapeA = x(2); % must be positive
% gammaBshapeB = x(3); % must be positive
p_doubleGamma = [6 8 1 1 0]; %guess
lb_DoubleGamma = [0 0 0 0 0];
ub_DoubleGamma = [15 15 15 1 10];
opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_doubleGamma, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFDoubleGamma, p_doubleGamma, t, hrf_Deconv, lb_DoubleGamma, ub_DoubleGamma, opts);
subplot(3,2,3)
% disp(x_doubleGamma)
fittedHRF = get_HRFDoubleGamma(x_doubleGamma,t);
plot(t,fittedHRF)
hold on
plot(t,get_HRFDoubleGamma(p_doubleGamma,t), '--r')
% deconvHRF = a(7,:)
deconvHRF = hrf_Deconv;
plot(t,deconvHRF)
title('Double Gamma')
legend('fitted params','starting params','deconv')

%%HRF GAMMA
% time = xdata;
% timelag = x(1);
% tau = x(2);
% exponent = x(3);
% amplitude = x(4);
% offset = x(5);
p_Gamma = [1 0.6 4 1 0]; %guess - Plot?
lb_Gamma = [0.1 0.1 2 0.1 0.1];
ub_Gamma = [16 inf 16 1 10];

opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_Gamma, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFGamma, p_Gamma, t, hrf_Deconv, lb_Gamma, ub_Gamma, opts);
subplot(3,2,4)
% disp(x_Gamma)
fittedHRF = get_HRFGamma(x_Gamma,t);
plot(t,fittedHRF)
hold on
plot(t,get_HRFGamma(p_Gamma,t), '--r')
% deconvHRF = a(7,:)
deconvHRF = hrf_Deconv;
plot(t,deconvHRF)
title('Gamma')
legend('fitted params','starting params','deconv')

% timelag = x(1);
% tau = x(2);
% exponent = x(3);
% timelag2 = x(4);
% tau2 = x(5);
% exponent2 = x(6);
% amplitude2 = x(7);
x_dGamma = [];
p_dGamma = [1 0.6 4 2 1.2 11 0.25]; %guess - Plot?
lb_dGamma = [0 0 0 0 0 0 -1];
ub_dGamma = [16 inf 16 16 inf 16 1];
opts = optimset('MaxFunEvals', 500, 'Display', 'off');
[x_dGamma, resnorm, ~, exitflag, output] = lsqcurvefit(@get_HRFDiffOfGamma, p_dGamma, t, hrf_Deconv, lb_dGamma, ub_dGamma, opts);
subplot(3,2,5)
% disp(x_dGamma)
fittedHRF = get_HRFDiffOfGamma(x_dGamma,t);
plot(t,fittedHRF)
hold on
plot(t,get_HRFDiffOfGamma(p_dGamma,t), '--r')
plot(t,hrf_Deconv)
title('Difference of Gamma')
legend('fitted params','starting params','deconv')

delayS = 2.5;
durationS = 2.5;
sampleDuration = 1.5;
totalDurationS = delayS+durationS;
totalDuration = round(totalDurationS/sampleDuration);  %total duration in samples
delay = round(delayS/sampleDuration);  %duration of delay in samples
duration = totalDuration - delay;
totalT = length(t);
hrfboxcar = [zeros(1,delay),ones(1,duration),zeros(1,totalT-totalDuration)];
subplot(3,2,6)
plot(t,hrfboxcar)
hold on
plot(t,hrf_Deconv)

title('Boxcar')
legend('Model','deconv')


%% Using Gramm to plot deconv HRF
% TR = 1.5;
% timePoints = 0 : TR : (TR * analysisParams.nHrfComponents - TR);

grammTimePoints = repmat(timePoints,size(HRF_TW_est,1),1);
grammFreqBin = repmat(freqBin',1,size(HRF_TW_est,2));

figure
g = gramm('x',grammTimePoints,'y',grammFreqBin,'z',HRF_TW_est);
g.geom_line();
% g.facet_grid([],grammFreqBin);
g.draw()

% figure
% clear g
% g = gramm('x',grammTimePoints,'y',HRF_TW_est)
% g.geom_line();
% % g.facet_wrap(grammFreqBin,[]);
% g.draw()

figure
clear g
g = gramm('x',grammTimePoints,'y',HRF_TW_est,'color',freqBin');
g.geom_line();
g.draw()

figure
clear g
g = gramm('x',grammTimePoints,'y',HRF_TW_est,'color',freqBin');
g.geom_line();
g.facet_grid(freqBin,[]);
g.draw()

end