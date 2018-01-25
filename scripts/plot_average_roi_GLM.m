function plot_average_roi_GLM(data,analName,stimInfo,plotInfo,figureName)
% plotLOGIC = what to plot based on number of stim

% eval(['roi_av = data.' analName '.roi_av;']);
% eval(['roi_pTW = data.' analName '.roi_pTW;']);
% eval(['roi_av_ratio = data.' analName '.roi_av_ratio;']);

if length(data.roi_av{1}) == stimInfo.sizes(1)
    xlabel = stimInfo.stimNames.bin;
    stimulusLevels = stimInfo.stimLevel_SL_bin;
    plotLOGIC = plotInfo.plotLOGIC.ROI_bin;
    nCols = 2;
elseif length(data.roi_av{1}) == stimInfo.sizes(2)
    xlabel = stimInfo.stimNames.mv;
    stimulusLevels = stimInfo.stimLevel_SL_mv;
    plotLOGIC = plotInfo.plotLOGIC.ROI_mv;
    nCols = 4;
else 
    xlabel = stimInfo.stimNames.all;
    stimulusLevels = stimInfo.stimLevel_SL;
    plotLOGIC = plotInfo.plotLOGIC.ROI_all;
    nCols = 4;
end

if plotLOGIC(1) == 1
%% compare ROI average Beta weights
plot_compareConditions_ROIav(data.roi_av{1},data.roi_av{2},xlabel,figureName)

end

% if plotLOGIC(2) == 1
% %% Compare average voxel tuning
% plot_compareConditions_pTW(roi_pTW{1}, roi_pTW{2},nCols,xlabel,figureName);
% 
% end

if plotLOGIC(3) == 1
    %% plot dB senstation level vs BOLD fMRI activity 
% [ fit ] =  plot_dbSLvsBetaWeight(ratio2Plot, level2Plot, fit , error2plot)
plot_dbSLvsBetaWeight(data.ratio2Plot,data.level2Plot,data.fit,data.error2plot,stimulusLevels);
end
end