function plot_splitAverage_roi_GLM(data,analName,stimInfo,plotInfo)
% plotLOGIC = what to plot based on number of stim

eval(['roi_av = data.' analName '.roi_av;']);
eval(['roi_pTW = data.' analName '.roi_pTW;']);

if length(roi_av{1}) == 8
    xlabel = stimInfo.stimNames.bin;
    plotLOGIC = plotInfo.plotLOGIC.ROI_bin;
    nCols = 2;
elseif length(roi_av{1}) == 28
    xlabel = stimInfo.stimNames.mv;
    plotLOGIC = plotInfo.plotLOGIC.ROI_mv;
    nCols = 4;
else 
    xlabel = stimInfo.stimNames.all;
    plotLOGIC = plotInfo.plotLOGIC.ROI_all;
    nCols = 4;
end

if plotLOGIC(1) == 1
%% compare ROI average Beta weights
plot_compareConditions_ROIav(roi_av{1},roi_av{2},xlabel)

end

if plotLOGIC(2) == 1
%% Compare average voxel tuning
plot_compareConditions_pTW(roi_pTW{1}, roi_pTW{2},nCols,xlabel);

end

end