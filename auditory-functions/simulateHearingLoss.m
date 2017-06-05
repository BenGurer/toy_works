function simulateHearingLoss
%
%   usage: simulateHearingLoss
%      by: Ben Gurer
%    date: 23/01/2017
% purpose: simulate hearing loss using masking noise
% 
% discription:
% Generate masking noise to simulate hearing loss. This function is intended to design
% a filter to be used in tdtMRI to generate masking noise which simulates hearing loss.

% Critical Ratio values taken from Hawkins And Stevens (1950) defined as:
% "Ratio between the monaural masked threshold of a pure tone and the level per cycle of the masking noise measured at the frequency of the pure tone".

% dBHL values taken from BS EN ISO 389-2:1997 (125Hz to 8kHZ) and ISO 389-5:2006 (8kHz to 16kHz).
% ERBs are calculated as per Glasberg and Moore (1990).
sampleDuration = 1/24.4140625;
f = 10.^(linspace(log10(20),log10(20000),100));

%% Create masking noise
% The masking noise aims to equalise the thresholds of hearing as a function of ERB. This means that stimuli presented above this masking noise will be of
% the same dB Sensation Level ie same level above threshold of perception.




%% Create hearing loss simulation masking noise
% This masking noise aims to simualte hearing loss by raising the threshold of auditory perception.

Threshold_dBHL = createSteeplySlopingHearingLoss_dBHL(f./1000);

RETSPL_int = getdBHLinSPL_inserts(f);

Threshold_dBSPL = Threshold_dBHL + RETSPL_int;

figure
plot(f,Threshold_dBSPL)
hold on
plot(f,Threshold_dBHL)
legend('Threshold in dB HL',...
    'Threshold in dB SPL',...
    'Location','best')
set(gca,'XLim',[min(f) max(f)])
title('Steeply sloping hearig loss according to screening crtiterion')
xlabel('Frequency (kHz)') 
ylabel('dB (HL/SPL)')

figure
plot(f,CriticalRatio_int,'k--')
hold on
plot(f,CriticalRatioPerERB,'r--')
plot(f,CriticalRatio_Zwicker,'b--')
legend('critical ratio (per Hz)',...
    'critical ratio (per ERB)',...
    'critical ratio (accroding to Zwicker)',...
    'Location','best')
xlabel('Frequency (kHz)') 
ylabel('dB')

figure
plot(f,Nee_lev,'k--')
hold on
plot(f,Nten_lev,'r--')
plot(f,Nee_CR_Zwicker_lev,'b--')

legend('EE filter',...
    'TEN filter (EE - CR)',...
    'EE filter - Zwicker CR',...
    'Location','best')
xlabel('Frequency (kHz)') 
ylabel('dB')


    
    lev = -10*log10(lcfErb(frq));
%     lev = lev/sqrt(mean(lev.^2));
    
    ee_1kHz = lev - (interp1(frq,lev,1));
    
    fftTotal = ee_1kHz + levelFFT;
    IntensityTotal = sum(10.^(fftTotal/10))*2;
    levelTotaldB = 10*log10(IntensityTotal);
    
    levelFFT = [levelFFT fliplr(levelFFT)];
    
% 
% figure
% plot(f,ERB)
