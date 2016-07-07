function d = toy_pRF
%(stimfile,TR)

% load one voxel time series - from auditory cortex with high R2

% desMatrix: stimulus presented at each time point

% create correlation matrix to find starting parameters
%     generate 10 by 10 matrix (f by TW) of possible pRF properties

thisView = getMLRView;
stimfile = viewGet(thisView,'stimfile',1);

% slice = 10;
% VoxelofInterest = (32,75,1,:);
% Load tSeries
tSeries = loadTSeries(thisView, 1, 'all');
VoxeltSeries = tSeries(26,72,13,:);
VoxeltSeries = squeeze(VoxeltSeries);

TR = 2;
stimNames = stimfile{1,1}.stimNames;
[nrows, ncols] = size(stimNames);

for k = 1:ncols
    StimulusSet(:,k) = sscanf(stimNames{:,k}, '%*s %d%*s', [1, inf]);
end;

StimulusSet = (StimulusSet)/1000;
nStimuli = length(stimfile{1,1}.stimNames);
tIndex = cell2mat(stimfile{1,1}.mylog.stimtimes_s);
tTime = length(VoxeltSeries);
dRaw = zeros(tTime,nStimuli);
desMatrix = zeros(tTime,nStimuli);
hrf = makeHrf(TR);


for i = 1:nStimuli
    dRaw(tIndex(i)+1:tIndex(i)+2,i) = 1; %+1;% add one indexing to handle stimulus at time 0
    dRaw(tIndex(i)+1,i) = 1;
    desMatrix(:,i) = conv(dRaw(:,i), hrf', 'same'); % same lenght as longest input)
end

figure
imagesc(desMatrix)

figure
plot(VoxeltSeries)

StimLowFreq = min(StimulusSet);
StimHighFreq = max(StimulusSet);
initalParams.pCF = lcfInvNErb(linspace(lcfNErb(StimLowFreq), lcfNErb(StimHighFreq), 10));
initalParams.pTW = [0.2 0.5 1 2];

% product of design matrix and pRF is the modelled time course
% corse search with correlation matrix
% use as initial parameters of fminsearch or whatever
modelT = makeModelledTimeCourse(desMatrix,StimulusSet,initalParams.pCF,initalParams.pTW);

r2 = zeros(length(initalParams.pCF),length(initalParams.pTW));

for i = 1:length(initalParams.pCF)
    for n = 1:length(initalParams.pTW)
        lm = fitlm(VoxeltSeries,squeeze(modelT(:,i,n)));
        r2(i,n) = lm.Rsquared.Ordinary;
    end
end

[minV, minI] = max(r2);
[minVV, minII] = max(minV);

INTpCF = initalParams.pCF(minI(minII));
INTpTW = initalParams.pTW (minII);

x0 = [INTpCF INTpTW 1];
xdata.desMatrix = desMatrix;
xdata.StimulusSet = StimulusSet;
ydata = double(VoxeltSeries);

modelData = pRFtimeCourse (x0,xdata);

lb = [0.02, 0.01];
ub= [20, 100];
fun = @pRFtimeCourse; % needs to be forward model
options = optimoptions('lsqcurvefit','TolFun',1e-10, 'MaxFunEvals', 400, 'MaxIter', 400); % default=1e-6
[x, error] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);

modelData = pRFtimeCourse (x,xdata);

figure
plot(modelData)
hold on
plot(ydata)

function hrf = makeHrf(TR)
% given the TR, return the HRF shape for t = 0 ... 30s
%
% using the equation given in the lecture (simple boynton version)

tau = 2; % time constant (s)
delta = 2; % time shift  (s)
t = [0:TR:30]; % vector of time points (in steps of TR)
% t = [0:1:30]; % vector of time points (in steps of seconds)

tshift = max(t-delta,0); % shifted, but not < 0
% hrf = (tshift ./ tau) .^2 .* exp(-tshift ./tau )./(2*tau);

hrf = thisGamma(t,1,delta,0,tau,3);

function gammafun = thisGamma(time,amplitude,timelag,offset,tau,exponent)

exponent = round(exponent);
% gamma function
gammafun = (((time-timelag)/tau).^(exponent-1).*exp(-(time-timelag)/tau))./(tau*factorial(exponent-1));

% negative values of time are set to zero,
% so that the function always starts at zero
gammafun(find((time-timelag) < 0)) = 0;

% normalize the amplitude
if (max(gammafun)-min(gammafun))~=0
  gammafun = (gammafun-min(gammafun)) ./ (max(gammafun)-min(gammafun));
end
gammafun = (amplitude*gammafun+offset);

function modelData = pRFtimeCourse (x,xdata)

pCF = x(1);
pTW = x(2);
scaling = x(3);

desMatrix = xdata.desMatrix;
StimulusSet = xdata.StimulusSet;

[nrows, ncols] = size(desMatrix);
modelData = zeros(nrows,1);

r = makepRF(pCF,pTW,StimulusSet);

% replace below with matrix *
% for k = 1:nrows
% modelData(k) = sum(desMatrix(k,:) .* r,2);
% end
modelData = desMatrix * r';
modelData = modelData./mean(modelData);
modelData = modelData .*scaling;

function modelTimeCourse = makeModelledTimeCourse(desMatrix,StimulusSet,pCF,pTW)

figure;
hold on
[nrows, ncols] = size(desMatrix);
% pRF = zeros(length(pCF),length(pTW),nrows);
modelTimeCourse = zeros(nrows,length(pCF),length(pTW));
for i = 1:length(pCF)
    % output vector length of stimulus set
    % use last arguement to set spacing of function
    for n = 1:length(pTW)
        r = makepRF(pCF(i),pTW(n),StimulusSet);
        % take pRF and multiply by design matrix
        % loop to lineraly sum
%         for k = 1:nrows
%             pRF(i,n,k) = sum(desMatrix(k,:) .* r,2);
%         end
        modelTimeCourse(:,i,n) = desMatrix * r';
        plot(modelTimeCourse(:,i,n))
    end
end
modelTimeCourse = modelTimeCourse./repmat(mean(modelTimeCourse),[nrows,1,1]);
    % change program to have time course on first dimention
% pRF = permute(pRF,[2,3,1]);
% take average of EACH time course - divide each time point by this average
% add 3rd param for fitting fuction - for scaling

function r = makepRF(pCF,pTW,stimf)
ERB = pTW*lcfErb(pCF);
P = 4*pCF/ERB;
g = abs(stimf-pCF)/pCF;
r = (1+P*g).*exp(-P*g);

function [w] = lcfROEX(param,f)
CF = param(1);
C =param(2);
P = 4*CF/C;
% ERB = C*lcfErb(CF);
% P = 4*CF/ERB;
g = abs(f-CF)/CF;
w = (1+P*g).*exp(-P*g);

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

