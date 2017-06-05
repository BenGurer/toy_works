
function checktdtOutputLevel

HB7Gain = -18;            % attenuation setting of the TDT HB7 Headphone Driver (in dB)
NNLGain = [-40 -27.4 -18.6 -12.21 -6.2 0]; %attenuation corresponding to the 6 NNL amplifier 'Acoustic Level' settings
NNLsetting = 6;              % amplification setting of NNL amplifier at the 7T scanner (in dB)
HB7CalibGain = -27;       % attenuation setting of the TDT HB7 Headphone Driver at which calibration was done
NNLCalibSetting = 6;         % amplification setting of NNL amplifier at the 7T scanner at which calibration was done

calibrationLevelLeft = 80.4; % calibration 04/03/2016 left side
calibrationLevelRight = 78.8; % calibration 04/03/2016 right side

% calibrationGainLeft level (in dB SPL) of a 1Volt 1kHz sinewave recorded at the inserts
calibrationGainLeft = - calibrationLevelLeft - 3; %corresponding level for a 1voltRMS noise
calibrationGainRight = - calibrationLevelRight - 3; %corresponding level for a 1voltRMS noise

calibrationGainLeft = calibrationGainLeft+HB7CalibGain-HB7Gain;

'NAmpL',10^((NLevel+calibrationGainLeft-LEE)/20)/scaling); %set the noise level

NLevel=eval(get(handleCaller,'String'))-SNR1dB;

NLevel = 35;%35;              % intended background noise level expressed in "masking" level (if Nlevel=Slevel, the signal would be just detectable)
%   SNR1dB = 10*log10(10^(1/10)-1); % SNR1dB is the SNR of a just detectable signal embedded in noise:
% a signal is just detectable if S+N is more intense than N alone by about 1 dB (Zwicker)
% solve 10*log10((IS+IN)/IN) = 1 dB for IS/IN and then apply 10*log10; IS/IN = signal/noise intensity;
% (i.e. a signal is just detectable in a noise level that's about 5.9 dB louder)
% this will be subtracted to the desired noise level (instead of adding it to the signal)
SNR1dB = 0;
%For the background noise, there is a twist because the intended sound level concerns the portion of the spectrum
%that stimulates one auditory (cochlear) filter and not the sound level of the total noise (at all frequencies)
%To account for this, we compute the ratio of the energy within one critical band around 1kHz and the total energy for an
%equally exciting noise (a noise that stimulates an auditory filter with the same energy)
LEE = lcfLEE(2^18,1,sampleDuration);
end


% ***** lcfLEE *****
function LEE = lcfLEE(N,F,sampleDuration)
%LEE is the level of an equally exciting stimulus with RMS of 0 dB (rms
%amplitude equals 1) with an ERB around F
DF = 1/(sampleDuration*N);
frq = DF*(1:N/2);

lev = -10*log10(lcfErb(frq));   %these two lines are equivalent to:
eeFilter = 10.^(lev/20);        % eeFilter = 1./sqrt(lcfErb(frq))   (see also lcfMakeNoise)

NF = lcfNErb(F);
F1 = lcfInvNErb(NF-0.5);
F2 = lcfInvNErb(NF+0.5);
bpFilter = zeros(1,N/2);
bpFilter(round(F1/DF):round(F2/DF)) = 1;

LEE = 10*log10(sum((eeFilter.*bpFilter).^2)/sum(eeFilter.^2));
end