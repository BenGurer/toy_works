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
