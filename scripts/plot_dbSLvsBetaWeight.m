function [ fit ] =  plot_dbSLvsBetaWeight(ratio2Plot, level2Plot, fit , error2plot,stimSensationLevel)

% ratio2Plot = [roi_av_ratio(stimulusLevels~=baseLevel_dB); mean(roi_av_ratio(stimulusLevels==baseLevel_dB))];
% 
% level2Plot = [stimulusLevels(stimulusLevels~=baseLevel_dB) baseLevel_dB];
% % 
% figure; scatter(level2Plot,ratio2Plot)
% figure; scatter(level2Plot_bin,ratio2Plot_bin)

%% dB SL vs beta weight
figure('color',[1 1 1]); 
scatter( level2Plot ,ratio2Plot )
hold on
% fit = polyfit(level2Plot',ratio2Plot,1);
plot(level2Plot,polyval(fit,level2Plot));
% correlation = corrcoef([level2Plot_mv' ratio2Plot_mv]);
% text(25,0.8,sprintf('Correlation = %.2f \n y = %.2fx + %.2f',correlation(2), fit(1), fit(2)))
% plot(f,f,'k--')
errorbar(error2plot(1),error2plot(2),error2plot(3))
xlabel('Stimulus Sensation Level (dB SL)'); ylabel('Average Beta Weight Ratio (ConB / ConA)');
plot(stimSensationLevel,stimSensationLevel./max(stimSensationLevel),'k--')

% figure; plot(stimInfo.stimLevel_SL_mv,stimInfo.stimLevel_SL_mv./max(stimInfo.stimLevel_SL_mv))

% 
% 
% figure; scatter(stimLevel_SL,peakWeightRatio)
% figure; scatter(stimLevel_SL_bin,peakWeightRatio_bin)
% figure; scatter(stimLevel_SL_mv,peakWeightRatio_mv)
% 
% figure; scatter(stimLevel_SL_bin,peakWeightDiff_bin)
% figure; scatter(stimLevel_SL_mv,peakWeightDiff_mv)

end