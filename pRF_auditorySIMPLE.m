function view = pRF_auditory(view,params)
%
% view = pRF_auditory(view,[params])
% 
% Loops throughs scans, loads corresponding tSeries, and computes time
% series statistics for each voxel:
% - mean 
% - median 
% - standard deviation
% - max frame-to-frame difference
% - max difference from median
%
% params: Optional initial parameters. Default: user is prompted via
%    GUI. Params must be a a structure with all of the following fields. 
% params.groupName: group of scans that will be analyzed.
%    Default: current group of view.
% params.scanList: vector specifying which scans to compute.
%    Default: all of the scans.
%
%
% Examples:
%
% paramss.groupName = 'Raw';
% n = viewGet([],'nScans',1)
% params.scanList = [1:n];
% view = pRF_auditory(view,params);
%
% view = pRF_auditory(view);
%
%
% djh, 7/2007
% $Id$	

if ~isview(view)
    help pRF_auditory
    mrErrorDlg('(pRF_auditory) Invalid view.')
end

% Get analysis parameters from pRF_auditoryGUI
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  groupName = viewGet(view,'groupName');
  params.groupName = groupName;
  % find out all the scans in the current group that have more than 2 volumes
  % since we can't compute tSeriesStats for scans with less
  groupNum = viewGet(view,'groupNum',groupName);
  for iScan = 1:viewGet(view,'nScans',groupNum)
    nVols(iScan) = viewGet(view,'nVolumes',iScan);
    if nVols(iScan) <= 2
      disp(sprintf('(pRF_auditory) Scan %s:%i has %i volumes - need at least 2 volumes to run tSeriesStats',groupName,iScan,nVols(iScan)));
    end
  end
  scanList = find(nVols>2);
  if length(scanList) > 0
    params.scanList = selectInList(view,'scans','Select scans for pRF_auditory',scanList);
    if isempty(params.scanList),return,end
  else
    disp(sprintf('(pRF_auditory) Could not find any scans to process'));
    return
  end
end

% Reconcile params with current status of group and ensure that params
% has the required fields.
params = defaultReconcileParams(params.groupName,params);

% Abort if params empty
if ieNotDefined('params')
  mrMsgBox('pRF_auditory cancelled',1);
  return
end

% Change group
groupName = params.groupName;
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',groupName);
if (groupNum ~= curGroup)
	mrWarnDlg(['Changing view to group: ',groupName]);
	view = viewSet(view,'currentGroup',groupNum);
end

% Compute it
[pRF_CF, pRF_TW] = computepRF_auditory(view,params);

% Make analysis structure
tsStats.name = 'pRF_auditory';  % This can be reset by editAnalysisGUI
tsStats.type = 'pRF_auditory';
tsStats.groupName = params.groupName;
tsStats.function = 'pRF_auditory';
tsStats.guiFunction = 'pRF_auditoryGUI';
tsStats.params = params;

% Install it in the view
view = viewSet(view,'newanalysis',tsStats);
view = viewSet(view,'newoverlay',pRF_CF);
view = viewSet(view,'newoverlay',pRF_TW);
view = viewSet(view,'newoverlay',tsStd);
view = viewSet(view,'newoverlay',tsMaxFrameDiff);
view = viewSet(view,'newoverlay',tsMaxMedianDiff);
view = viewSet(view,'newoverlay',tsMeanDividedByStd);

% Save it
saveAnalysis(view,tsStats.name);


keyboard


return

function [pRF_CF, pRF_TW] = computepRF_auditory(view,params)

% Get nScans from view and get scanList from params
scanList = params.scanList;
nScans = viewGet(view,'nScans');

% Intialize overlay structures

% CF
pRF_CF.name = 'CF';
pRF_CF.function = 'pRF_auditory';
pRF_CF.data = cell(1,nScans);
pRF_CF.params = params;
pRF_CF.colormap = jet(256);
pRF_CF.groupName = params.groupName;
pRF_CF.interrogator = 'pRF_auditoryPlot';

% TW
pRF_TW.name = 'TW';
pRF_TW.function = 'pRF_auditory';
pRF_TW.data = cell(1,nScans);
pRF_TW.params = params;
pRF_TW.colormap = jet(256);
pRF_TW.groupName = params.groupName;
pRF_TW.interrogator = 'pRF_auditoryPlot';

disp('Computing pRF...');
warning('off','MATLAB:divideByZero');
for scanIndex=1:length(scanList)
   scanNum = scanList(scanIndex);
   waitHandle = mrWaitBar(0,['Computing pRF for Scan ' int2str(scanNum) ':']);

   % sliceDims: [ydim xdim] for single slice
   % volDims; [ydim xdim nslices] for single scan
   sliceDims = viewGet(view,'sliceDims',scanNum);
   volDims = viewGet(view,'dims',scanNum);

   % Initialize data with NaNs
   pRF_CF.data{scanNum} = NaN*ones(volDims);
   pRF_TW.data{scanNum} = NaN*ones(volDims);


   nslices = viewGet(view,'nslices',scanNum);
     for sliceNum = 1:nslices
         % modify to loop through each voxel here or in side comptepRF
     [pRF_CFSeries,pRF_TWSeries] = computepRF_auditorySeries(view,scanNum,sliceNum);
     pRF_CF.data{scanNum}(:,:,sliceNum) = reshape(pRF_CFSeries,sliceDims);
     pRF_TW.data{scanNum}(:,:,sliceNum) = reshape(pRF_TWSeries,sliceDims);

     % Update waitbar
     mrWaitBar(sliceNum/nslices,waitHandle);
   end
   mrCloseDlg(waitHandle);
end

% Fill range fields 
pRF_CF.range = findRange(pRF_CF.data);
pRF_TW.range = findRange(pRF_TW.data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pRF_CFSeries,pRF_TWSeries] = ...
  computepRF_auditorySeries(view,scan,slice)

% Get junk frames and nframes
junkframes = viewGet(view,'junkframes',scan);
nframes = viewGet(view,'nframes',scan);

% Load tSeries
tSeries = loadTSeries(view, scan, slice);

% Reshape the tSeries
tSeries = reshapeTSeries(tSeries);

% Remove junkFrames
tSeries = tSeries(junkframes+1:junkframes+nframes,:);

%% This is where we calculated the pRF

% tsMeanSeries = nanmean(tSeries);
% tsMedianSeries = nanmedian(tSeries);
% tsStdSeries = nanstd(tSeries);
% tsMaxFrameDiffSeries = nanmax(abs(tSeries(2:end,:)-tSeries(1:end-1,:)));
% tsMaxMedianDiffSeries = nanmax(abs(tSeries - repmat(tsMedianSeries,[nframes,1])));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range = findRange(data)

ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
     thisData = data{scan}(:);
     thisData = thisData(~isinf(thisData));
    ampMin = min([ampMin min(thisData)]);
    ampMax = max([ampMax max(thisData)]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end
