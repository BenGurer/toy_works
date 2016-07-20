function [pCF,pTW,error] = pRFpush (VoxeltSeries,stiminfo,TR,stimTR)

if any(isnan(VoxeltSeries))
    
    pCF = nan;
    pTW = nan;
    error = nan;
else

% make data one row double
VoxeltSeries = squeeze(VoxeltSeries);
VoxeltSeries = percentTSeries(VoxeltSeries,'detrend','Linear','spatialNormalization','Divide by mean','subtractMean', 'Yes', 'temporalNormalization', 'No');
 


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
coarseSearch = makeModelledTimeCourse(stiminfo.desMatrix,stiminfo.StimulusSet,initalParams.pCF,initalParams.pTW);

% Calculate correlation of each modelled time course with actual voxel time series
r2 = zeros(length(initalParams.pCF),length(initalParams.pTW));
for i = 1:length(initalParams.pCF)
    for n = 1:length(initalParams.pTW)
        lm = fitlm(VoxeltSeries,squeeze(coarseSearch.modelTimeSeries(:,i,n)));
        r2(i,n) = lm.Rsquared.Ordinary;
    end
end

% find pRF model parameters with highest correlation and uses are inital
% parameters for minimising search
r2(r2~=r2)=0;
[maxV, maxI] = max(r2);
[maxVV, maxII] = max(maxV);
INTpCF = initalParams.pCF(maxI(maxII));
INTpTW = initalParams.pTW (maxII);

% set minimising search parameters
x0 = [INTpCF INTpTW 1];
xdata.desMatrix = stiminfo.desMatrix;
xdata.StimulusSet = stiminfo.StimulusSet;
ydata = double(VoxeltSeries);
lb = [0.02, 0.0001];
ub= [20, 1000];
fun = @pRFtimeCourse; % needs to be forward model
[x, error] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);

pCF = x(1);
pTW = x(2);
end
% % Modelled time series using highest correlated parameters from coarse search
% initialpRFmodelTimeSeries = pRFtimeCourse (x0,xdata);
% % Modelled time series using minimising search parameters
% fittedpRFmodelTimeSeries = pRFtimeCourse (x,xdata);

% f = [0.02:1:20];
% r = makepRF(x(1),x(2),f);

%%%%%%%%%%%%%%%%%%
%% Plot figures %%
%%%%%%%%%%%%%%%%%%

% figure; plot(f,r)
% 
% %% Hemodynamic Response Function
% figure; plot(hrf); title('HDR function'), xlabel('Time (TR)'), ylabel('Amplitude')
% 
% %% Design matrix
% figure
% imagesc(desMatrix); title('Design Matrix'), xlabel('Stimulus ID'), ylabel('Time (TR)')
% 
% %% Voxel Time Series
% % figure
% % plot(VoxeltSeries); title('Voxel Time Series'), xlabel('Time (TR)'), ylabel('BOLD Percent Signal Change')
% 
% %% Modelled time courses
% % Transpose to allow waterfall plot but otherwise not used
% a = permute(coarseSearch.modelTimeSeries,[2,1,3]);
% [szx, szy, szz] = size(a);
% figure;
% for i = 1:szz
%     subplot(round(szz/2),round(szz/2),i)
%     waterfall(a(:,:,i));title('Model Time Series'), xlabel('Time (TR)'), ylabel('Initial Charactersic Frequency ID'), zlabel('Amplitude')
% end
% 
% %% Modelled pRF
% figure;
% for i = 1:length(initalParams.pCF)
%     subplot(round(length(initalParams.pCF)/2),round(length(initalParams.pCF)/2),i)
%     for n = 1:length(initalParams.pTW)
%         plot(coarseSearch.pRF(:,i,n))
%         hold on
%     end
% end
% title('Modelled pRFs'), xlabel('Frequency (ERB Scale)'), ylabel('Amplitude')
% %% R2 values of fitted modelled time series
% figure; contourf(r2); colorbar; title('R2 values of Model fits'), xlabel('Initial Tuning Width ID'), ylabel('Initial Charactersic Frequency ID')
% 
% %% Voxel time series and fitted model time series
% figure
% hold on
% plot(ydata)
% plot(initialpRFmodelTimeSeries)
% plot(fittedpRFmodelTimeSeries); title('Voxel and Modelled Time Series'), xlabel('Time (TR)'), ylabel('BOLD Percent Signal Change')
% text(max(xlim)-range(xlim)/10,max(ylim)-range(ylim)/20,sprintf('pCF = %.2f kHz\npTW = %.2f ERB',x(1),x(2)))
% keyboard

% hrf = hrf./ max(hrf);

function modelTimeCourse = pRFtimeCourse (x,xdata)
pCF = x(1);
pTW = x(2);
scaling = x(3);
% offset = x(3);

desMatrix = xdata.desMatrix;
StimulusSet = xdata.StimulusSet;

[nrows, ncols] = size(desMatrix);
% modelTimeCourse = zeros(nrows,1);

% Response of voxel to each stimulus
r = makepRF(pCF,pTW,StimulusSet);

% Matrix multiplication of design matrix and response of voxel to each
% stimulus =
% At each time point(colums) sum the product of each element in the row (stimulus
% presenation convolved with HRF) with the voxels pRF response to each stimulus
modelTimeCourse = desMatrix * r';
% modelTimeCourse = modelTimeCourse + offset; % off set to account for data off set
modelTimeCourse = modelTimeCourse./mean(modelTimeCourse); % give data a mean of 1
modelTimeCourse = modelTimeCourse .* scaling; % scale of data unknown so account for this


function coarseSearch = makeModelledTimeCourse(desMatrix,StimulusSet,pCF,pTW)
[nrows, ncols] = size(desMatrix);
% pRF = zeros(length(pCF),length(pTW),nrows);
coarseSearch.modelTimeSeries = zeros(nrows,length(pCF),length(pTW));
coarseSearch.pRF  = zeros(ncols,length(pCF),length(pTW));
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
        coarseSearch.modelTimeSeries(:,i,n) = desMatrix * r';
        %         plot(modelTimeCourse(:,i,n))
        coarseSearch.pRF(:,i,n)=r;
    end
end
coarseSearch.modelTimeSeries = coarseSearch.modelTimeSeries +1;
coarseSearch.modelTimeSeries = coarseSearch.modelTimeSeries./repmat(mean(coarseSearch.modelTimeSeries),[nrows,1,1]); % take average of EACH time course - divide each time point by this average


function r = makepRF(pCF,pTW,stimf)
% take into account bandwidth of noise
% ERB = pTW*lcfErb(pCF);
% P = 4*pCF/ERB;
P = 4*pCF/pTW;
g = abs(stimf-pCF)/pCF;
r = (1+P*g).*exp(-P*g);

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

