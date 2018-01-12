function data = overlay_difference(overlayA,overlayB)

r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'pRF_auditory';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(v,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(v,'nScans'));
r2.range = [0 1];
r2.clip = [0 1];
% colormap is made with a little bit less on the dark end
r2.colormap = hot(312);
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'pRFPlot_HDRauditory';
r2.mergeFunction = 'pRFMergeParams';

    r2.data{scanNum} = nan(scanDims);
    
            r2.data{scanNum}(x(i),y(i),z(i)) = thisr2(i);
             r2.params{scanNum} = thisParams;
 
  thisView = viewSet(thisView,'newoverlay',outputOverlay);            
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