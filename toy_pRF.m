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
VoxeltSeries = tSeries(26,72,13,:); % r2 = 0.311 max index = 7 tw index 7ish
% VoxeltSeries = tSeries(32,72,10,:); % r2 = 0.312 maxIndex 9 TW 4 - 5
VoxeltSeries = squeeze(VoxeltSeries);

TR = 2;
stimNames = stimfile{1,1}.stimNames;
[nrows, ncols] = size(stimNames);

for k = 1:ncols
    StimulusSet(:,k) = sscanf(stimNames{:,k}, '%*s %d%*s', [1, inf]);
end;

% need to divid by TR to get time points in TR
% change for loop with +1 used before
% get stim times for second run
StimulusSet = (StimulusSet)/1000;
nStimuli = length(stimfile{1,1}.stimNames);
% tIndex = cell2mat(stimfile{1,1}.mylog.stimtimes_s);
tTime = length(VoxeltSeries);
tIndex = [cell2mat(stimfile{1,1}.mylog.stimtimes_s)/2;(cell2mat(stimfile{1,2}.mylog.stimtimes_s)/2)+tTime/2];
tIndex = tIndex+1; % add one indexing to handle stimulus at time 0
dRaw = zeros(tTime,nStimuli);
desMatrix = zeros(tTime,nStimuli);
hrf = makeHrf(TR);


for i = 1:nStimuli
%     dRaw(tIndex(i)+1:tIndex(i)+2,i) = 1; %+1;% add one indexing to handle stimulus at time 0
%     dRaw(tIndex(i)+1,i) = 1;

    dRaw(tIndex(1,i),i) = 1;
    dRaw(tIndex(2,i),i) = 1;
    
    desMatrix(:,i) = conv(dRaw(:,i), hrf', 'same'); % same lenght as longest input)
end

figure
imagesc(desMatrix); title('Design Matrix'), xlabel('Stimulus ID'), ylabel('Time (TR)')

figure
plot(VoxeltSeries); title('Voxel Time Series'), xlabel('Time (TR)'), ylabel('BOLD Percent Signal Change')

StimLowFreq = min(StimulusSet);
StimHighFreq = max(StimulusSet);
initalParams.pCF = lcfInvNErb(linspace(lcfNErb(StimLowFreq), lcfNErb(StimHighFreq), 10));
initalParams.pTW = [0.5 1 5 10 50 100];

% product of design matrix and pRF is the modelled time course
% corse search with correlation matrix
% use as initial parameters of fminsearch or whatever
modelT = makeModelledTimeCourse(desMatrix,StimulusSet,initalParams.pCF,initalParams.pTW);
a = permute(modelT,[2,1,3]);
[szx szy szz] = size(a);
figure;
for i = 1:szz
    subplot(round(szz/2),round(szz/2),i)
 waterfall(a(:,:,i));title('Model Time Series'), xlabel('Time (TR)'), ylabel('Initial Charactersic Frequency ID'), zlabel('Amplitude')
% figure; waterfall(modelT(:,:,1))
end

r2 = zeros(length(initalParams.pCF),length(initalParams.pTW));

for i = 1:length(initalParams.pCF)
    for n = 1:length(initalParams.pTW)
        lm = fitlm(VoxeltSeries,squeeze(modelT(:,i,n)));
        r2(i,n) = lm.Rsquared.Ordinary;
    end
end

r2(r2~=r2)=0;
[maxV, maxI] = max(r2);
[maxVV, maxII] = max(maxV);


INTpCF = initalParams.pCF(maxI(maxII));
INTpTW = initalParams.pTW (maxII);

figure; contourf(r2); colorbar; title('R2 values of Model fits'), xlabel('Initial Tuning Width ID'), ylabel('Initial Charactersic Frequency ID')

x0 = [INTpCF INTpTW 1 1];
xdata.desMatrix = desMatrix;
xdata.StimulusSet = StimulusSet;
ydata = double(VoxeltSeries);

modelData = pRFtimeCourse (x0,xdata);

lb = [0.02, 0.0001];
ub= [20, 1000];
fun = @pRFtimeCourse; % needs to be forward model
% options = optimoptions('lsqcurvefit','TolFun',1e-10, 'MaxFunEvals', 400, 'MaxIter', 400); % default=1e-6
% [x, error] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
[x, error] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);

modelData = pRFtimeCourse (x,xdata);

figure
hold on
plot(ydata)
plot(modelData); title('Voxel and Modelled Time Series'), xlabel('Time (TR)'), ylabel('BOLD Percent Signal Change')
text(max(xlim)-range(xlim)/10,max(ylim)-range(ylim)/20,sprintf('pCF = %.2f kHz\npTW = %.2f ERB',x(1),x(2)))

function hrf = makeHrf(TR)
% given the TR, return the HRF shape for t = 0 ... 30s
%
t = [0:TR:30]; % vector of time points (in steps of TR)
x = 4;
y = 11;
z = 4;
hrf = gampdf(t,x,1)-gampdf(t,y,1)/z;
figure; plot(hrf); title('HDR function'), xlabel('Time (TR)'), ylabel('Amplitude')

% function hrf = makeHrf(TR)
% using the equation given in the lecture (simple boynton version)
% 
% tau =1.5; % time constant (s)
% delta = 3; % time shift  (s)
% t = [0:TR:30]; % vector of time points (in steps of TR)
% % t = [0:1:30]; % vector of time points (in steps of seconds)
% 
% tshift = max(t-delta,0); % shifted, but not < 0
% % hrf = (tshift ./ tau) .^2 .* exp(-tshift ./tau )./(2*tau);
% 
% hrf = thisGamma(t,1,delta,0,tau,3);

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
offset = x(4);

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
modelData = modelData + offset;
modelData = modelData./mean(modelData);
modelData = modelData .* scaling;


function modelTimeCourse = makeModelledTimeCourse(desMatrix,StimulusSet,pCF,pTW)



[nrows, ncols] = size(desMatrix);
% pRF = zeros(length(pCF),length(pTW),nrows);
modelTimeCourse = zeros(nrows,length(pCF),length(pTW));
pRF  = zeros(ncols,length(pCF),length(pTW));
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
        %         plot(modelTimeCourse(:,i,n))
        pRF(:,i,n)=r;
    end
end
% modelTimeCourse = modelTimeCourse +1;
% modelTimeCourse = modelTimeCourse./repmat(mean(modelTimeCourse),[nrows,1,1]);

figure;
for i = 1:length(pCF)
    for n = 1:length(pTW)
        plot(modelTimeCourse(:,i,n))
        hold on
    end
end

figure;
for i = 1:length(pCF)
    subplot(round(length(pCF)/2),round(length(pCF)/2),i)
    for n = 1:length(pTW)
        plot(pRF(:,i,n))
        hold on
    end
end

% title('Modelled Time Courses'), xlabel('Time (TR)'), ylabel('Amplitude')

    % change program to have time course on first dimention
% pRF = permute(pRF,[2,3,1]);
% take average of EACH time course - divide each time point by this average
% add 3rd param for fitting fuction - for scaling

function r = makepRF(pCF,pTW,stimf)
% take into account bandwidth of noise
% ERB = pTW*lcfErb(pCF);
% P = 4*pCF/ERB;
P = 4*pCF/pTW;
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

