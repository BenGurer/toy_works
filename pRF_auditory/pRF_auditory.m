% pRF_auditory.m
%
%        $Id:$ 
%      usage: pRF_auditory(v,params,varargin)
%         by: only slightly modified from pRF.m by justin gardner
%       date: 
%    purpose: compute pRF analysis on MLR data
%
%             if you just want a default parameter structure you
%             can do:
% 
%             v = newView;
%             [v params] = pRF_auditory(v,[],'justGetParams=1','defaultParams=1','scanList=1')
%
%             Note that justGetParams,defualtParams and scanList are independent parameters, so
%             if you want, say to bring up the GUI to set the params, but not run the analysis, you
%             can do:
%             [v params] = pRF_auditory(v,[],'justGetParams=1');
%


% TO DO
% Add HDR estimates to corse search
% Create function to plot average HDR estimates of ROI
% give option to fit in ERB space
% make HDR fitting a tick box in GUI - to change case in pRF_auditoryFit
% error check -
%       - does it still work with continuous?
%       - check pRF_auditoryFit to do list
% space initial params on model scale
% calculate the sigma search range. min = response to one stimulus. max =
% respond to all equally/broadband
% change to stimulus domain - just change when creating stim.x
% could add option to set hrf lendth in seconds
function [v d] = pRF_auditory(v,params,varargin)

% check arguments
if nargin < 1
  help pRF_auditory
  return
end

d = [];
% a version number in case we make major changes
pRFVersion = 1;

% params defaults to empty
if nargin < 2,params =[];end

% other arguments
justGetParams=[];defaultParams=[];scanList=[];
groupNum=[];
getArgs(varargin,{'justGetParams=0','defaultParams=0','scanList=[]','groupNum=[]'});

% first get parameters
if isempty(params)
  % get group
  if isempty(groupNum),groupNum = viewGet(v,'curGroup');end
  % put up the gui
  params = pRFGUI_auditory('v',v,'groupNum',groupNum,'defaultParams',defaultParams,'scanList',scanList);
end

% just return parameters
if justGetParams,d = params;return,end

% Reconcile params with current status of group and ensure that it has
% the required fields. 
params = defaultReconcileParams([],params);

% Abort if params empty
if isempty(params),return,end

% check the params
params = checkPRFparams(params);

% set the group
v = viewSet(v,'curGroup',params.groupName);

% create the parameters for the r2 overlay
dateString = datestr(now);
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

% at this point we need to decide on which parameters we want to estimate
% from data
%
% for pRF_auditory e.g.
%       PrefCentreFreq (1, 2, 3)
%       rfHalfWidth...?
%       etc.

% create the parameters for the PrefCentreFreq overlay
PrefCentreFreq = r2;
PrefCentreFreq.name = 'PrefCentreFreq';
PrefCentreFreq.range = [0.02 20]; %cf 0.5-3.5
PrefCentreFreq.clip = [0.02 20];
PrefCentreFreq.colormapType = 'normal';
PrefCentreFreq.colormap = jet(256);

if any(strcmp(params.pRFFit.voxelScale,{'lin'}))
elseif any(strcmp(params.pRFFit.voxelScale,{'log'}))
pCFscaled = r2;
pCFscaled.name = 'pCFscaled';
pCFscaled.range = [log10(0.02) log10(20)];
pCFscaled.clip = [log10(0.02) log10(20)];
pCFscaled.colormapType = 'normal';
pCFscaled.colormap = jet(256);
elseif any(strcmp(params.pRFFit.voxelScale,{'erb'}))
pCFscaled = r2;
pCFscaled.name = 'pCFscaled';
pCFscaled.range = [funNErb(0.02) funNErb(20)];
pCFscaled.clip = [funNErb(0.02) funNErb(20)];
pCFscaled.colormapType = 'normal';
pCFscaled.colormap = jet(256);
else
  disp(sprintf('(pRFFit) Unknown voxelScale: %s',fitParams.voxelScale));
end


% create the paramteres for the rfHalfWidth overlay
% deal with the sigma.
rfHalfWidth = r2;
rfHalfWidth.name = 'rfHalfWidth';
rfHalfWidth.range = [0 100];
rfHalfWidth.clip = [0 100];
rfHalfWidth.colormapType = 'setRangeToMax';
rfHalfWidth.colormap = jet(256);

% create the parameters for the HDR exponent overlay
hdrExp = r2;
hdrExp.name = 'hdrExponent';
hdrExp.range = [0 15];
hdrExp.clip = [0 16];
hdrExp.colormapType = 'setRangeToMax';
hdrExp.colormap = jet(256);

% create the parameters for the NRMSD overlay
NRMSD = r2;
NRMSD.name = 'NRMSD';
NRMSD.range = [0 1];
NRMSD.clip = [-30 30];
NRMSD.colormapType = 'setRangeToMax';
NRMSD.colormap = hot(312);
% colormap is made with a little bit less on the dark end
NRMSD.colormap = NRMSD.colormap(end-255:end,:);

% create the paramteres for the hdrtimelag overlay
hdrtimelag = r2;
hdrtimelag.name = 'hdrtimelag';
hdrtimelag.range = [0 15];
hdrtimelag.clip = [0 16];
hdrtimelag.colormapType = 'setRangeToMax';
hdrtimelag.colormap = hot(256);

% create the paramteres for the hdrScale overlay
hdrScale = r2;
hdrScale.name = 'hdrScale';
hdrScale.range = [0 100];
hdrScale.clip = [0 100];
hdrScale.colormapType = 'setRangeToMax';
hdrScale.colormap = jet(256);

% get number of workers 
nProcessors = mlrNumWorkers;

% code snippet for clearing precomputed prefit
%global gpRFFitStimImage;gpRFFitStimImage = [];

dispHeader
disp(sprintf('(pRF_auditory) Running on scans %s:%s (restrict %s)',params.groupName,num2str(params.scanNum,'%i '),params.restrict ));

for scanNum = params.scanNum
  % see how long it took
  tic;
  
  % get voxels that we are restricted to
  [x y z] = getVoxelRestriction(v,params,scanNum);
  if isempty(x)
    disp(sprintf('(pRF_auditory) No voxels to analyze with current restriction'));
    return
  end

  % get total number of voxels
  n = length(x);

  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum);
  
%   if params.pRFFit.supersampling == 1   
      var.supersamplingMode = 'Automatic';
      params.pRFFit.d = getStimvolpRF(v,var);
%   end
  
  % init overlays
  r2.data{scanNum} = nan(scanDims);
  PrefCentreFreq.data{scanNum} = nan(scanDims);
  pCFscaled.data{scanNum} = nan(scanDims);
  NRMSD.data{scanNum} = nan(scanDims);
  rfHalfWidth.data{scanNum} = nan(scanDims);
  hdrExp.data{scanNum} = nan(scanDims);
  hdrtimelag.data{scanNum} = nan(scanDims);
  hdrScale.data{scanNum} = nan(scanDims);

  % default all variables that will be returned
  % by pRFFIt, so that we can call it the
  % second time and save some time
  concatInfo = [];
  stim = [];
  
  % save pRF parameters
  pRFAnal.d{scanNum}.ver = pRFVersion;
  pRFAnal.d{scanNum}.linearCoords = [];
  pRFAnal.d{scanNum}.params = [];
  pRFAnal.d{scanNum}.fitParams = [];
  
  % get some information from pRFFit that will be used again in
  % the fits, including concatInfo, stim, prefit, etc.
  fit = pRF_auditoryFit(v,scanNum,[],[],[],'fitTypeParams',params.pRFFit,'returnPrefit',true);
  if isempty(fit),return,end
  stim = fit.stim;
  pRFAnal.d{scanNum}.stim = cellArray(stim);
  pRFAnal.d{scanNum}.stimX = fit.stimX;
  pRFAnal.d{scanNum}.stimY = fit.stimY;
  pRFAnal.d{scanNum}.stimT = fit.stimT;
  concatInfo = fit.concatInfo;
  pRFAnal.d{scanNum}.concatInfo = fit.concatInfo;
  prefit = fit.prefit;
  paramsInfo = fit.paramsInfo;
  pRFAnal.d{scanNum}.paramsInfo = paramsInfo;
  pRFAnal.d{scanNum}.fitParams = fit.fitParams;
  % grab all these fields and stick them onto a structure called paramsInfo
  % preallocate some space
  rawParams = nan(fit.nParams,n);
  r = nan(n,fit.concatInfo.n);
  thisr2 = nan(1,n);
  thisNRMSD = nan(1,n);
  thisRfHalfWidth = nan(1,n);
  thishdrExp = nan(1,n);
  thishdrtimelag = nan(1,n);
  thishdrscale = nan(1,n);

  % get some info about the scan to pass in (which prevents
  % pRFFit from calling viewGet - which is problematic for distributed computing
  framePeriod = viewGet(v,'framePeriod');
  junkFrames = viewGet(v,'junkFrames',scanNum);

  % compute pRF for each voxel in the restriction
  if params.pRFFit.prefitOnly,algorithm='prefit-only';else,algorithm=params.pRFFit.algorithm;end

  % disp info about fitting
  dispHeader;
  disp(sprintf('(pRF_auditory) Scan %s:%i (restrict %s) running on %i processor(s)',params.groupName,scanNum,params.restrict,nProcessors));
  disp(sprintf('(pRF_auditory) Computing %s fits using %s for %i voxels',params.pRFFit.rfType,algorithm,n));
  dispHeader;

  % this is a bit arbitrary but is the number of voxels to read in at a time.
  % should probably be either calculated based on memory demands or a
  % user settings. The bigger the number the less overhead and will run faster
  % but consume more memory. The overhead is not terribly significant though
  % as tested on my machine - maybe a few percent faster with full n, but
  % on many machines without enough memory that will crash it so keeping
  % this preliminary value in for now.
  blockSize = n;
  tic;
  % break into blocks of voxels to go easy on memory
  % if blockSize = n then this just does on block at a time.
  for blockStart = 1:blockSize:n

    % display information about what we are doing
    % get blockEnd
    blockEnd = min(blockStart + blockSize-1,n);
    blockSize = blockEnd-blockStart+1;
    
    % load ROI
    loadROI = makeEmptyROI(v,'scanNum',scanNum,'groupNum',params.groupName);
    loadROI.coords(1,1:blockSize) = x(blockStart:blockEnd);
    loadROI.coords(2,1:blockSize) = y(blockStart:blockEnd);
    loadROI.coords(3,1:blockSize) = z(blockStart:blockEnd);
    % load all time series for block, we do this to pass into pRFFit. Generally
    % the purpose here is that if we run on distributed computing, we
    % can't load each voxel's time series one at a time. If this is
    % too large for memory then you can comment this out and not
    % pass it into pRFFit and pRFFit will load the tSeries itself
    loadROI = loadROITSeries(v,loadROI,scanNum,params.groupName);
    %     blockEnd = loadROI.n;
    %     blockSize = loadROI.n;
    n = loadROI.n;
    blockEnd = min(blockStart + blockSize-1,n);
    blockSize = blockEnd-blockStart+1;
    % reorder x,y,z coordinates since they can get scrambled in loadROITSeries
    x(blockStart:blockEnd) = loadROI.scanCoords(1,1:blockSize);
    y(blockStart:blockEnd) = loadROI.scanCoords(2,1:blockSize);
    z(blockStart:blockEnd) = loadROI.scanCoords(3,1:blockSize);
    % keep the linear coords
    pRFAnal.d{scanNum}.linearCoords = [pRFAnal.d{scanNum}.linearCoords sub2ind(scanDims,x(blockStart:blockEnd),y(blockStart:blockEnd),z(blockStart:blockEnd))];
    
    if blockStart ~= 1
        % display time update
      dispHeader(sprintf('(pRF_auditory) %0.1f%% done in %s (Estimated time remaining: %s)',100*blockStart/n,mlrDispElapsedTime(toc),mlrDispElapsedTime((toc*n/blockStart) - toc)));
    end

    % now loop over each voxel
%     for i = blockStart:blockEnd
    parfor i = blockStart:blockEnd
      fit = pRF_auditoryFit(v,scanNum,x(i),y(i),z(i),'stim',stim,'concatInfo',concatInfo,'prefit',prefit,'fitTypeParams',params.pRFFit,'dispIndex',i,'dispN',n,'tSeries',loadROI.tSeries(i-blockStart+1,:)','framePeriod',framePeriod,'junkFrames',junkFrames,'paramsInfo',paramsInfo);
      if ~isempty(fit)
        % keep data, note that we are keeping temporarily in
        % a vector here so that parfor won't complain
        % then afterwords we put it into the actual overlay struct
        thisr2(i) = fit.r2;
        thisPrefCentreFreq(i) = fit.PrefCentreFreq;
        thispCFscaled(i) = fit.pCFscaled;
        thisRfHalfWidth(i) = fit.rfHalfWidth;
        % keep parameters
        rawParams(:,i) = fit.params(:);
        thishdrExp(i) = fit.hdrExp;
        thishdrtimelag(i) = fit.hdrtimelag;
        thishdrscale(i) = fit.scale(1);
        thishdroffset(i) = fit.scale(2);
        r(i,:) = fit.r;
        thisNRMSD(i) = fit.NRMSD;
        
      end
    end
      
    % set overlays
    for i = 1:n
      r2.data{scanNum}(x(i),y(i),z(i)) = thisr2(i);
      NRMSD.data{scanNum}(x(i),y(i),z(i)) = thisNRMSD(i);
      PrefCentreFreq.data{scanNum}(x(i),y(i),z(i)) = thisPrefCentreFreq(i);
      pCFscaled.data{scanNum}(x(i),y(i),z(i)) = thispCFscaled(i);
      rfHalfWidth.data{scanNum}(x(i),y(i),z(i)) = thisRfHalfWidth(i);
      hdrExp.data{scanNum}(x(i),y(i),z(i)) = thishdrExp(i);
      hdrtimelag.data{scanNum}(x(i),y(i),z(i)) = thishdrtimelag(i);
      hdrScale.data{scanNum}(x(i),y(i),z(i)) = thishdrscale(i);
    end
  end
  % display time update
  dispHeader;
  disp(sprintf('(pRF_auditory) Fitting %i voxels took %s.',n,mlrDispElapsedTime(toc)));
  dispHeader;
  
  pRFAnal.d{scanNum}.params = rawParams;
  pRFAnal.d{scanNum}.r = r;
  pRFAnal.d{scanNum}.r2 = thisr2;
  
  iScan = find(params.scanNum == scanNum);
  thisParams.scanNum = params.scanNum(iScan);
  r2.params{scanNum} = thisParams;
  NRMSD.params{scanNum} = thisParams;
  PrefCentreFreq.params{scanNum} = thisParams;
  pCFscaled.params{scanNum} = thisParams;
  rfHalfWidth.params{scanNum} = thisParams; 
  hdrExp.params{scanNum} = thisParams; 
  hdrtimelag.params{scanNum} = thisParams; 
  hdrScale.params{scanNum} = thisParams; 
  
  % display how long it took
  disp(sprintf('(pRF_auditory) Fitting for %s:%i took in total: %s',params.groupName,scanNum,mlrDispElapsedTime(toc)));
  
end


% install analysis
pRFAnal.name = params.saveName;
pRFAnal.type = 'pRFAnal';
pRFAnal.groupName = params.groupName;
pRFAnal.function = 'pRF_auditory';
pRFAnal.reconcileFunction = 'defaultReconcileParams';
pRFAnal.mergeFunction = 'pRFMergeParams';
pRFAnal.guiFunction = 'pRF_auditoryGUI';
pRFAnal.params = params;
pRFAnal.overlays = [r2 NRMSD PrefCentreFreq pCFscaled rfHalfWidth hdrExp hdrtimelag hdrScale];
pRFAnal.curOverlay = 1;
pRFAnal.date = dateString;
v = viewSet(v,'newAnalysis',pRFAnal);

% if we are going to merge, temporarily set overwritePolicy
if isfield(params,'mergeAnalysis') && params.mergeAnalysis
  saveMethod = mrGetPref('overwritePolicy');
  mrSetPref('overwritePolicy','Merge');
end
% Save it
saveAnalysis(v,pRFAnal.name);
% now set policy back
if isfield(params,'mergeAnalysis') && params.mergeAnalysis
  mrSetPref('overwritePolicy',saveMethod);
end

if ~isempty(viewGet(v,'fignum'))
  refreshMLRDisplay(viewGet(v,'viewNum'));
end

%set(viewGet(v,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(d)
    pRFAnal.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(pRFAnal.d) == 1
    d = pRFAnal.d{1};
  else
    d = pRFAnal.d;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getVoxelRestriction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x y z] = getVoxelRestriction(v,params,scanNum)

x = [];y = [];z = [];

if strncmp(params.restrict,'Base: ',6)
  % get the base name
  baseName = params.restrict(7:end);
  baseNums = [];
  if strcmp(baseName,'ALL')
    for iBase = 1:viewGet(v,'numBase')
      % if the base is a surface or flat then add to the list
      if any(viewGet(v,'baseType',iBase) == [1 2])
	baseNums(end+1) = iBase;
      end
    end
  else
    baseNums = viewGet(v,'baseNum',baseName);
  end
  % cycle through all bases that we are going to run on
  scanCoords = [];
  for iBase = 1:length(baseNums)
    % get the baseNum
    baseNum = baseNums(iBase);
    if isempty(baseNum)
      disp(sprintf('(pRF_auditory) Could not find base to restrict to: %s',params.restrict));
      continue
    end
    % get the base
    base = viewGet(v,'base',baseNum);
    if isempty(base)
      disp(sprintf('(pRF_auditory) Could not find base to restrict to: %s',params.restrict));
      return;
    end
    % if flat or surface
    if any(base.type == [1 2])
      % get base coordinates from the coordMap
      for corticalDepth = 0:0.1:1
	if base.type == 1
	  % flat map
	  baseCoords = (base.coordMap.innerCoords + corticalDepth * (base.coordMap.outerCoords-base.coordMap.innerCoords));
	  baseCoords = reshape(baseCoords,prod(size(base.data)),3)';
	else
	  % surface
	  baseCoords = (base.coordMap.innerVtcs + corticalDepth * (base.coordMap.outerVtcs-base.coordMap.innerVtcs))';
	end
	% convert to 4xn array
	baseCoords(4,:) = 1;
	% and convert to scan coordinates
	base2scan = viewGet(v,'base2scan',scanNum,params.groupName,baseNum);
	scanCoords = [scanCoords round(base2scan*baseCoords)];
      end
    end
  end
  % check against scandims
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  scanCoords = mrSub2ind(scanDims,scanCoords(1,:),scanCoords(2,:),scanCoords(3,:));
  % remove duplicates and nans
  scanCoords = scanCoords(~isnan(scanCoords));
  scanCoords = unique(scanCoords);
  % convert back to x,y,z coordinates
  [x y z] = ind2sub(scanDims,scanCoords);
elseif strncmp(params.restrict,'ROI: ',5)
  % get the roi name
  roiName = params.restrict(6:end);
  scanCoords = getROICoordinates(v,roiName,scanNum,params.groupName,'straightXform=1');
  if isempty(scanCoords),return,end
  x = scanCoords(1,:);y = scanCoords(2,:);z = scanCoords(3,:);
elseif strncmp(params.restrict,'None',4)
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  [x y z]  = ndgrid(1:scanDims(1),1:scanDims(2),1:scanDims(3));
  x = x(:);y = y(:);z = z(:);
else
  keyboard
end

%check if we have already computed Voxels
if isfield(params,'computedVoxels') && (length(params.computedVoxels)>=scanNum) && ~isempty(params.computedVoxels{scanNum})
  % get scan dims
  scanDims = viewGet(v,'scanDims',scanNum,params.groupName);
  % convert x, y, z to linear coords
  linearCoords = sub2ind(scanDims,x,y,z);
  % get new ones
  newLinearCoords = setdiff(linearCoords,params.computedVoxels{scanNum});
  if length(newLinearCoords) ~= length(linearCoords)
    % show what we are doing
    disp(sprintf('(pRF) Dropping %i voxels that have been already computed',length(linearCoords)-length(newLinearCoords)));
    % convert back to x, y, z
    [x y z] = ind2sub(scanDims,newLinearCoords);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%    checkPRFparams    %
%%%%%%%%%%%%%%%%%%%%%%%%
function params = checkPRFparams(params)


% check the pRFFit params
checkFields = {{'stimImageDiffTolerance',5}};
for iFit = 1:length(params.pRFFit)

  % set defaults
  for iField = 1:length(checkFields)
    if ~isfield(params.pRFFit(iFit),checkFields{iField}{1})
      params.pRFFit(iFit).(checkFields{iField}{1}) = checkFields{iField}{2};
    end
  end
end

