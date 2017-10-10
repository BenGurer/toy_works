function stimData = sHL_setupStudyParams
%% get stimulus properties
%% get stimulus frequencies
lowFreqkHz = 0.1;
highFreqkHz = 8;
nStim = 32;
[stimData.stimFreqs, stimData.stimFreqs_bin, stimData.stimFreqs_mv] = convertStimIDtoFrequency(lowFreqkHz,highFreqkHz,nStim);

%% get stimulus senssation level
[stimData.stimLevel_SL, maskingLevel] = calStimulusSensationLevel(stimFreqs);

figure('color',[1 1 1]); plot(stimFreqs,stimLevel_SL);xlabel('Frequency (kHz)'); ylabel('Sensation Level (dB SL)')

% bin sensation level
stimLevel_SL_bin = [];
binSize = 4;
c = 1;
for i = 1:length(stimLevel_SL)/binSize
    stimLevel_SL_bin(i) = mean(stimLevel_SL(c:c + (binSize-1)));
    c = c + binSize;
end
stimData.stimLevel_SL_bin = stimLevel_SL_bin;


% moving average sensation level
nBins = 8;
windowAvSize = length(stimLevel_SL)/nBins;

if isreal(windowAvSize) && rem(windowAvSize,1)==0
    
    loopLength = length(stimLevel_SL) - windowAvSize;
    stimLevel_SL_mv = zeros(1,loopLength);
    
    for i = 1:loopLength
        stimLevel_SL_mv(i) = mean(stimLevel_SL(i:i+windowAvSize-1));
    end
    
else
    error('Moving average window not an integer')
end

stimData.stimLevel_SL_mv = stimLevel_SL_mv;