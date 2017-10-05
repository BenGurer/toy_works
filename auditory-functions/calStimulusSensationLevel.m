function [stimLevel_SL, maskingLevel] = calStimulusSensationLevel(stimFreqs)

threshold = 75;
stimLevel_SPL = 75;
baselineLevel = 25;
hearingLossSim = funSimulateHearingLoss(stimFreqs);
hearingLossSim = min(hearingLossSim,threshold);
Baseline = baselineLevel * ones(size(stimFreqs));
maskingLevel = max(hearingLossSim,Baseline);
stimLevel_SL = stimLevel_SPL - maskingLevel;