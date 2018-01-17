function [ thisView , differenceData ] = script_createDifferenceMaps(thisView,overlayA,overlayB)

differenceData = abs(overlayB.data{1} - overlayA.data{1});
scanDims = size(overlayB.data{1});
% create overlay structure
difference = overlayB;
dateString = datestr(now);
r2.date = dateString;
difference.name = 'difference';
% difference.colorRange = [0 1];
% difference.range = [0 1];
% difference.clip = [0 1];
difference.range = [min(min(min(differenceData))) max(max(max(differenceData)))];
difference.clip = [min(min(min(differenceData))) max(max(max(differenceData)))];
difference.colorRange = [min(min(min(differenceData))) max(max(max(differenceData)))];
if exist('brewermap.m', 'file')
    difference.colormap = brewermap(256,'Blues');
else
    difference.colormap = hot(256);
end
% difference.alpha = 1;
% difference.colormapType = 'setRangeToMax';
difference.data{1} = nan(scanDims);
difference.data{1} = differenceData;
% difference.data{scanNum}(x(i),y(i),z(i)) = thisr2(i);

thisView = viewSet(thisView,'newoverlay',difference,4);
%              'newoverlay'
%    view = viewSet(view,'newoverlay',overlayStructure,[analysisNum]);
%
%    val must be a structure or an array of structures with fields
%      - name: string
%      - groupName: string group name
%      - range: [min max] values, possibly excluding meaningless values
%      - params: structure specifying arguments to function
%          To recompute: view = function(view,params)
%    and optionnally
%      - function: string function name by which it was computed
%      - reconcileFunction: string function name that reconciles params
%         and data with tseries file names.
%         [newparams,newdata] = reconcileFunction(groupName,params,data)
%      - interrogator: function that gets called when you click on the
%          overlay (e.g., to produce a time series plot).
%      - data: cell array of [x y z] arrays
%      - date: specifies when it was computed, typically generated using
%          datestr(now)
%      - clip: [min max] to be displayed/thresholded
%      - colormap: 256x3 array of RGB values
%      - alpha: transparency value for alphaSlider
%      - colorRange: [min max] of the colormap


end