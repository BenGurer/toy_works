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

f = 10.^(linspace(log10(20),log10(20000),10));

%% Create masking noise
% The masking noise aims to equalise the thresholds of hearing as a function of ERB. This means that stimuli presented above this masking noise will be of
% the same dB Sensation Level ie same level above threshold of perception.

CriticalRatio_int = getCriticalRatio_HawkinsAndStevens1950(f);
ERB = getERB(f/1000);

CriticalRatioPerERB = CriticalRatio_int-10*log10(ERB);
CriticalRatio_Zwicker = 10*log10(10^(1/10)-1)*ones(size(f));

Nee_lev = -10*log10(ERB);
Nten_lev = -10*log10(ERB)-CriticalRatioPerERB;
Nee_CR_Zwicker_lev = -10*log10(ERB)-CriticalRatio_Zwicker;


%% Create hearing loss simulation masking noise
% This masking noise aims to simualte hearing loss by raising the threshold of auditory perception.

Threshold_dBHL = createSteeplySlopingHearingLoss_dBHL(f./1000);

RETSPL_int = getdBHLinSPL_inserts(f);

Threshold_dBSPL = Threshold_dBHL + RETSPL_int;

figure
plot(f,Threshold_dBSPL)
hold on
plot(f,RETSPL_int)

figure
plot(f,CriticalRatio_int,'r')
hold on
plot(f,CriticalRatioPerERB,'b')
plot(f,CriticalRatio_Zwicker,'g')

figure
plot(f,Nee_lev,'b')
hold on
plot(f,Nten_lev,'r')
plot(f,Nee_CR_Zwicker_lev,'g')

% 
% figure
% plot(f,ERB)
