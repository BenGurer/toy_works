function thisView = pRF_auditory_oneVoxel
%(view,params)
%
% view = pRF_auditory(view,[params])
%
% Loops throughs scans, loads corresponding tSeries, and population receptive field for each voxel:
% - population charactersic frequency (pCF)
% - population tuning width (pTW)
% - fit error (error)
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
% bjg, 8/2016
% $Id$

% load current mrLoadRet view
thisView = getMLRView;

if ~isview(thisView)
    help pRF_auditory
    mrErrorDlg('(pRF_auditory) Invalid view.')
end

% Get analysis parameters from pRF_auditoryGUI
if ieNotDefined('params')
    % Initialize analysis parameters with default values
    groupName = viewGet(thisView,'groupName');
    params.groupName = groupName;
    % find out all the scans in the current group that have more than 2 volumes
    % since we can't compute tSeriesStats for scans with less
    groupNum = viewGet(thisView,'groupNum',groupName);
    for iScan = 1:viewGet(thisView,'nScans',groupNum)
        nVols(iScan) = viewGet(thisView,'nVolumes',iScan);
        if nVols(iScan) <= 2
            disp(sprintf('(pRF_auditory) Scan %s:%i has %i volumes - need at least 2 volumes to run tSeriesStats',groupName,iScan,nVols(iScan)));
        end
    end
    scanList = find(nVols>2);
    if length(scanList) > 0
        params.scanList = selectInList(thisView,'scans','Select scans for pRF_auditory',scanList);
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
curGroup = viewGet(thisView,'currentGroup');
groupNum = viewGet(thisView,'groupNum',groupName);
if (groupNum ~= curGroup)
    mrWarnDlg(['Changing view to group: ',groupName]);
    thisView = viewSet(thisView,'currentGroup',groupNum);
end

% Compute it
[pRF_CF, pRF_TW] = computepRF_auditory(thisView,params);

% Make analysis structure
pRF.name = 'pRF_auditory';  % This can be reset by editAnalysisGUI
pRF.type = 'pRF_auditory';
pRF.groupName = params.groupName;
pRF.function = 'pRF_auditory';
pRF.guiFunction = 'pRF_auditoryGUI';
pRF.params = params;

% Install it in the view
thisView = viewSet(thisView,'newanalysis',pRF);
thisView = viewSet(thisView,'newoverlay',pRF_CF);
thisView = viewSet(thisView,'newoverlay',pRF_TW);

% Save it
saveAnalysis(thisView,pRF.name);

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
    
    %% Load Time Series
    % Get junk frames and nframes
    junkframes = viewGet(view,'junkframes',scanNum);
    nframes = viewGet(view,'nframes',scanNum);
    
    % Load tSeries
    tSeries = loadTSeries(view, scanNum, 'all');
    
    % Reshape the tSeries
    % tSeries = reshapeTSeries(tSeries);
    %     % Remove junkFrames
    %     tSeries = tSeries(junkframes+1:junkframes+nframes,:);
    % Total time of data set
    %     tTime = size(tSeries,4); % time in TR
    
    % Load stimulus file
    stimfile = viewGet(view,'stimfile',scanNum);
    
%     [d stimNames var] = getStimvol(view,stimVariable,varargin);
    
    stimTR = 2;
    TR = 2;
    
    stiminfo = makeStimInfo(stimfile,nframes,stimTR,TR);
    
    %%%%%%%%%%%%%%%%%%
    %% Coarse search %%
    %%%%%%%%%%%%%%%%%%
    % product of design matrix and pRF is the modelled time course
    % coarse search of pRF using correlation matrix
    % Best fit from this is used as initial parameters for minimising search
    
    % search between limits of stimulus set frequency range
    StimLowFreq = min(stiminfo.StimulusSet);
    StimHighFreq = max(stiminfo.StimulusSet);
    initalParams.pCF = lcfInvNErb(linspace(lcfNErb(StimLowFreq), lcfNErb(StimHighFreq), 10));
    initalParams.pTW = [0.5 1 5 10 50 100];
    
    % Returns a matrix of modelled time courses and pRF
    coarseSearch = makeModelledTimeCourse(stiminfo.designMatrix,stiminfo.StimulusSet,initalParams.pCF,initalParams.pTW);
    
    [pCF,pTW,error] = computepRF_auditorySeries(tSeries(25,67,12,:),stiminfo,coarseSearch,initalParams);
    %     nslices = viewGet(view,'nslices',scanNum);
    %     for sliceNum = 1:nslices
    %         [pCF,pTW,error] = computepRF_auditorySeries(tSeries(:,:,sliceNum,:),stiminfo,coarseSearch,initalParams);
    %         pRF_CF.data{scanNum}(:,:,sliceNum) = pCF;
    %         pRF_TW.data{scanNum}(:,:,sliceNum) = pTW;
    %         pRF_error.data{scanNum}(:,:,sliceNum) = error;
    %         % Update waitbar
    %         mrWaitBar(sliceNum/nslices,waitHandle);
    %     end
    mrCloseDlg(waitHandle);
end

% Fill range fields
pRF_CF.range = findRange(pRF_CF.data);
pRF_TW.range = findRange(pRF_TW.data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stiminfo = makeStimInfo(stimfile,nframes,stimTR,TR)
% stiminfo.StimulusSet
% stiminfo.designMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make Hemodynamic Response Function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrf = makeHrf(TR);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Design Matrix %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get stimulus information &
stimNames = stimfile{1,1}.stimNames; % Load Stimulus names
[nrows, ncols] = size(stimNames);
StimulusSet = zeros(nrows, ncols);
for k = 1:ncols
    StimulusSet(:,k) = sscanf(stimNames{:,k}, '%*s %d%*s', [1, inf]); % remove text to get frequency in Hz
end
stiminfo.StimulusSet = (StimulusSet)/1000; % convert frequency from Hz to kHz
nStimuli = length(stimfile{1,1}.stimNames); % Number of stimuli presented

tIndex = [cell2mat(stimfile{1,1}.mylog.stimtimes_s)/TR;(cell2mat(stimfile{1,2}.mylog.stimtimes_s)/TR)+(max(cell2mat(stimfile{1,1}.mylog.stimtimes_s))+2)/TR]; % Load stim times. Divide by TR to get time points in TR
tIndex = tIndex+1; % add one indexing to handle stimulus at time 0
dRaw = zeros(nframes,nStimuli); % Stimulus Impulse Matrix
stiminfo.designMatrix = zeros(nframes,nStimuli); % Design matrix - Stimulus impulse convolved with HRF
for i = 1:nStimuli
    dRaw(tIndex(1,i),i) = 1;
    dRaw(tIndex(2,i),i) = 1;
    stiminfo.designMatrix(:,i) = conv(dRaw(:,i), hrf', 'same'); % same lenght as longest input
end

function hrf = makeHrf(TR)
% given the TR, return the HRF shape for t = 0 ... 30s
%
% t = [0:TR:16./TR]; % vector of time points (in steps of TR)
t = [0:1:16./TR]; % vector of time points (in steps of TR)
x = 4/TR;
y = 11/TR;
z = 4;
hrf = gampdf(t,x,1)-gampdf(t,y,1)/z;

function [pCF,pTW,error] = computepRF_auditorySeries(tSeries,stiminfo,coarseSearch,initalParams)

%% This is where we calculated the pRF
% Calculate TR from number of stim file time points compared to image timepoints?
% loop through slice voxels
% output 3 matrixs of same size with pCF, pTW and error of each voxel
[nrows, ncols, nslice, nframes] = size(tSeries);
% nslice should be 1 - doing 1 slice at a time

pCF = zeros(nrows,ncols);
pTW = zeros(nrows,ncols);
error = zeros(nrows,ncols);
for i = 1:nrows
    for ii = 1:ncols
        VoxeltSeries = tSeries(i,ii,1,:);
        [pCF(i,ii),pTW(i,ii),error(i,ii)] = pRFpush (VoxeltSeries,stiminfo,coarseSearch,initalParams);
    end
end
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
function erb = lcfErb(f)
% ***** lcfErb *****
% ERBs as per Glasberg and Moore (1990);
A = 24.7/1000; B = 4.37;
erb = A*(B*f+1);

function nerb = lcfNErb(f)
% ***** lcfNErb *****
% Converts frequency to ERB number;
A = 24.7/1000; B = 4.37;
nerb = 1/(A*B)*log(B*f+1);

function f = lcfInvNErb(nerb)
% ***** lcfInvNErb *****
% Converts ERB number to frequency;
A = 24.7/1000; B = 4.37;
f = 1/B*(exp(A*B*nerb)-1);