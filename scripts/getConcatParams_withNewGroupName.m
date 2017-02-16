
function params = getConcatParams_withNewGroupName(view,newGroupName,varargin)

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end


interpTypes = {'nearest','linear','spline','cubic'};
if find(strcmp(mrGetPref('interpMethod'),interpTypes))
  interpTypes = putOnTopOfList(mrGetPref('interpMethod'),interpTypes);
end

% check to see if any of the scans in this group have a non identity scan2scan
needToWarp = 0;
for i = 2:viewGet(view,'nScans')
  if ~isequal(viewGet(view,'scan2scan',1,[],i),eye(4))
    needToWarp = 1;
  end
end

% check to see if any scans have a tSense that is not one, if it
% is then we will offer to notch filter the data, otherwise hide the option
offerNotchFilter = false;
defaultNotchFilterSetting = false;
if ~isempty(viewGet(view,'groupNum','Raw'))
  for iScan = 1:viewGet(view,'nScans','Raw')
    tSense = viewGet(view,'auxParam','tSense',iScan,'Raw');
    if iscell(tSense),tSense = cell2mat(tSense);end
    if isscalar(tSense) && (tSense > 1)
      offerNotchFilter = true;
      % set default setting to false if the current scan is a concat
      if ~isempty(viewGet(view,'concatInfo'))
        defaultNotchFilterSetting = false;
      end
      break;
    end
  end
end

% description of paramaters (used by mrParamsDialog functions)
paramsInfo = {...
    {'groupName',putOnTopOfList(viewGet(view,'groupName'),viewGet(view,'groupNames')),'Name of group from which to make concatenation'},...
    {'newGroupName',newGroupName,'Name of group to which concatenation will be saved. If group does not already exist, it will be created.'},...
    {'description','Concatenation of [x...x]','Description that will be set to have the scannumbers that are selected'},...
    {'filterType',{'Detrend and highpass','Detrend only','None'},'Which filter to use. Highpass filtering will use the cutoff below.'},...
    {'filterCutoff',0.01,'minmax=[0 inf]','Highpass filter cutoff in Hz'},...
    {'percentSignal',1,'type=checkbox','Convert to percent signal change. This is done by simply dividing by the mean, so that you get a timecourse where the mean is 1. (The mean is not subtracted out so that if subsequent analysis tries to divide by mean again it will not affect the time series)'}};
if needToWarp
  paramsInfo{end+1} = {'warp',1,'type=checkbox','Warp images based on alignment. This can be used to concatenate together scans taken on different days. If the scans have the same slice prescription this will not do anything.'};
  paramsInfo{end+1} = {'warpInterpMethod',interpTypes,'Interpolation method for warp','contingent=warp'};
end
paramsInfo{end+1} = {'projectOutMeanVector',0,'type=checkbox','Project out a mean vector defined from one roi out of the data. This is used if you want to remove a global component that is estimated as the mean over some roi from your data. If you select this you will be asked to choose one roi for defining the mean vector of what you want to project out and another roi from which you want to project out.'};
if offerNotchFilter
  paramsInfo{end+1} = {'notchFilterForTSense',defaultNotchFilterSetting,'type=checkbox','This is used to notch out the highest frequency for tSense data'};
end
    
% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    params = mrParamsDialog(paramsInfo);
  end
  % no params means user hit cancel
  if isempty(params),return,end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select scans
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  view = viewSet(view, 'groupName', params.groupName);
  if ~ieNotDefined('scanList')
    params.scanList = scanList;
  elseif defaultParams
    params.scanList = 1:viewGet(view,'nScans');
  else
    params.scanList = selectInList(view,'scans','Select Scans',1:viewGet(view,'nScans'));
  end
  if isempty(params.scanList),return,end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % if warp is set, then ask which scan to use as base scan for warp, unless set by default
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~needToWarp,params.warp = 0;params.warpInterpMethod = interpTypes{1};end
  if params.warp
    if defaultParams
      params.warpBaseScan = params.scanList(1);
    else
      % create a list of scans to ask user which scan to use as the base scan
      for i = 1:viewGet(view,'nScans')
	scanNames{i} = sprintf('%i:%s (%s)',i,viewGet(view,'description',i),viewGet(view,'tSeriesFile',i));
      end
      warpParams = mrParamsDialog({{'warpBaseScan',scanNames,'The scan that will be used as the base scan to warp all the other scans to'}});
      if isempty(warpParams),return,end
      params.warpBaseScan = find(strcmp(warpParams.warpBaseScan,scanNames));
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Get rois for doing projection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if params.projectOutMeanVector || defaultParams
    % projectOutMeanVectorParams
    params.projectOutMeanVectorParams = projectOutMeanVectorParams(view,defaultParams);
    if isempty(params.projectOutMeanVectorParams),return,end
  end
  % check the parameters
  params = mrParamsReconcile(params.groupName,params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params.paramsInfo = paramsInfo;
  params = mrParamsReconcile(params.groupName,params);
end