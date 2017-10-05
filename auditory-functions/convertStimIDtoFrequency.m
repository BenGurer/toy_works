function [stimFreqs, stimFreqs_bin, stimFreqs_mv] = convertStimIDtoFrequency(lowFreqkHz,highFreqkHz,nStim)

stimNERBs = funNErb(lowFreqkHz):(funNErb(highFreqkHz)-funNErb(lowFreqkHz))/(nStim-1):funNErb(highFreqkHz);

stimFreqs = funInvNErb(stimNERBs);

% bin stimulus

stimFreqs_bin = [];
binSize = 4;
c = 1;
for i = 1:length(stimNERBs)/binSize
stimFreqs_bin(i) = funInvNErb(mean(stimNERBs(c:c + (binSize-1))));
c = c + binSize;
end

% moving average

nBins = 8;
windowAvSize = length(stimNERBs)/nBins;

if isreal(windowAvSize) && rem(windowAvSize,1)==0
    
loopLength = length(stimNERBs) - windowAvSize;
stimFreqs_mv = zeros(1,loopLength);

for i = 1:loopLength
stimFreqs_mv(i) = funInvNErb(mean(stimNERBs(i:i+windowAvSize-1)));
end

else
    error('Moving average window not an integer')
end

