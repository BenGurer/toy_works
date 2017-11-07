function [ratio, level, fit, error] = cal_dbSLvsBetaWeight(roi_av_ratio,stimulusLevels,baseLevel_dB)

ratio = [roi_av_ratio(stimulusLevels~=baseLevel_dB); mean(roi_av_ratio(stimulusLevels==baseLevel_dB))];

level = [stimulusLevels(stimulusLevels~=baseLevel_dB) baseLevel_dB];

fit = polyfit(level',ratio,1);

error = [mean(stimulusLevels(stimulusLevels==50)), mean(roi_av_ratio(stimulusLevels==50)), std(roi_av_ratio(stimulusLevels==50))];

end