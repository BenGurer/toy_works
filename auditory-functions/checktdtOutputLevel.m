function checktdtOutputLevel

HB7Gain = -18;            % attenuation setting of the TDT HB7 Headphone Driver (in dB)
HB7CalibGain = -27;       % attenuation setting of the TDT HB7 Headphone Driver at which calibration was done

% S14 calibrationLevel = level in dBSPL when 1kHz sine wave presented at 1Volt
% calibrationGainLeft level (in dB SPL) of a 1Volt 1kHz sinewave recorded at the inserts
calibrationLevelLeft = 80.4; % calibration 04/03/2016 left side
calibrationLevelRight = 78.8; % calibration 04/03/2016 right side

% make negative because it is going to be used to normalised the output level ie set to 0dB
% - 3 to make 1Volt RMS rather than peak
calibrationGainLeft = - calibrationLevelLeft - 3; %corresponding level for a 1voltRMS noise
calibrationGainRight = - calibrationLevelRight - 3; %corresponding level for a 1voltRMS noise

% take into account the gain settings of the HB7 headphone driver - both when used for calibration and for when the signal is presented (now)
calibrationGainLeft = calibrationGainLeft+HB7CalibGain-HB7Gain;
     
% the desired level of a stimulus want the masking noise to mask.
NLevel = 35;%35;              % intended background noise level expressed in "masking" level (if Nlevel=Slevel, the signal would be just detectable)
%   SNR1dB = 10*log10(10^(1/10)-1); % SNR1dB is the SNR of a just detectable signal embedded in noise:
% a signal is just detectable if S+N is more intense than N alone by about 1 dB (Zwicker)
% solve 10*log10((IS+IN)/IN) = 1 dB for IS/IN and then apply 10*log10; IS/IN = signal/noise intensity;
% (i.e. a signal is just detectable in a noise level that's about 5.9 dB louder)
% this will be subtracted to the desired noise level (instead of adding it to the signal)

%For the background noise, there is a twist because the intended sound level concerns the portion of the spectrum
%that stimulates one auditory (cochlear) filter and not the sound level of the total noise (at all frequencies)
%To account for this, we compute the ratio of the energy within one critical band around 1kHz and the total energy for an
%equally exciting noise (a noise that stimulates an auditory filter with the same energy)
sampleDuration = 0.041
LEE = lcfLEE(2^18,1,sampleDuration); %-15.6923

maxVoltage = 10; %saturation voltage of TDT
scaling = maxVoltage*0.9/max(max(fNoise));
%write background noise to TDT
invoke(RP2,'WriteTagVEX','FNoise',0,'I16',round(scaling*fNoise/10*2^15));  % fill the noise buffer with 16-bit integers  


invoke(RP2,'SetTagVal','NAmpL',10^((calibrationGainLeft-LEE)/20)/scaling); %set the Nom noise level - the level is now set in the creation of the noise, lcfMakeNoise ----> lcfSetNoiseLevel
invoke(RP2,'SetTagVal','NAmpR',10^((calibrationGainRight-LEE)/20)/scaling); %set the Nom noise level
invoke(RP2,'SetTagVal','SplitScale',maxVoltage/(2^15-1)); %set the scaling factor that converts the signals from 16-bit integers to floats after splitting the two channels

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