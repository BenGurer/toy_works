function pRFdev_compressive

% Get stim weighting
% convert to intnesity
% compare

% create gaussian pRF
% apply power law
% plot both

x = 0:0.25:20;
mu = 4;
sigma = 0.5;
compressFun = 10;

stimLevel = 75;
nMaskingLevel = 25;
stimSLlevel = stimLevel - nMaskingLevel;
Threshold_dBHL = funSimulateHearingLoss(x);
thresholdBaseline = nMaskingLevel*ones(size(x));
threshEvel =  min(max(Threshold_dBHL,thresholdBaseline),stimLevel);
stimWeightingdBSL = stimLevel-threshEvel;

%     Convert to pressure
% Lp(dB SPL) = 20 log10 p/p0
% p0 = 0.00002 pa
% p(Pa) = p0 .10.^ Lp(dB SPL)/20
stimWeightingPressure = (2*10.^-5) .* (10.^(stimWeightingdBSL/20));

stimWeightingIntensity = (2*10.^-5) .* (10.^(stimWeightingdBSL/10));
stimWeightingdBSLCheck = 20 .* log10(stimWeightingPressure./(2*10.^-5));
% stimWeightingIntensity = 10.^-12 .* 10.^((20 .* log10((stimWeightingPressure/10.^-5))./10));
% stimWeightingIntensity = (stimWeightingPressure.^2) ./ 400;

% stimWeightingPressureCheck = sqrt(stimWeightingIntensity.*400);

stimWeightingPressureCheck = sqrt(stimWeightingIntensity);

stimWeightingSIL = 10 .* (log10((stimWeightingIntensity./(2*10.^-5))));
% stimWeightingSIL = 10 .* (log10((stimWeightingIntensity./(10.^-12))));
stimWeightingIntensityNorm = stimWeightingIntensity/max(stimWeightingIntensity);
stimWeightingdBSLNorm = stimWeightingdBSL/max(stimWeightingdBSL);

rfModel = exp(-(((x-mu).^2)/(2*(sigma^2))));

figure('Color','white');

subplot(2,2,4)
plot(x,rfModel,'b')
hold on
plot(x,rfModel.^compressFun,'r')
plot(x,max(1+20.*log10(rfModel),0),'g')
plot(x,(5.*log10(rfModel))/max(5.*log10(rfModel)),'g--')
title('pRF model')

subplot(2,2,1)
plot(x,stimWeightingdBSLNorm)
hold on
plot(x,stimWeightingdBSLNorm.^2)
plot(x,stimWeightingIntensityNorm)
plot(x,stimWeightingIntensityNorm.^.2)
title('Stimulus Weighting - Normalised')
legend('dB SL','dB SL^2','Intensity','Intensity ^2')

subplot(2,2,2)
plot(x,stimWeightingdBSL)
hold on
plot(x,stimWeightingSIL)
plot(x,stimWeightingdBSLCheck)
title('Stimulus Weighting dB')

subplot(2,2,3)
yyaxis left
plot(x,stimWeightingPressure)
hold on
yyaxis left
plot(x,stimWeightingPressureCheck)
ylabel('Pressure')

yyaxis right
plot(x,stimWeightingIntensity)
ylabel('Intensity')
title('Stimulus Weighting - Pressure and Intensity')
