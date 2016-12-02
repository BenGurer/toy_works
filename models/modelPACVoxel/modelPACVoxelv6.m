function [Ctrd,Sprd,dCtrd,dSprd,eCtrd,eSprd,CF,C,fitCF,fitC,fitTW,fiterror] = modelPACVoxelv2
%% THE QUESTION

% what question I am going to answer and what I will present in my report

% In x time how can we best estimate the properties of a voxel?
% experimental variables:
    % number of repeats
    % number of conditions
        % sweeps or discrete
    % TR time
        % Sparse = 480 scans per 60 mins
        % 3.5 = ~1029 scans per 60 mins
% experimental constants
    % total duration
    % presentable frequency range

% data set varibles:
   % CF
   % C
   % noise
   
% how many scans can we collect in experiement time
% how to use these scans - more:
% repeats
% conditions

%% TO DO
% convert C to be linear on different scales - change input to that domain then convert output back to frequency 
% print what each figure is, and graph
% if no maximum - error message or something

%% NOTES
% fwhm = 3;% spread of hemodynamic response in mm
% sigma = fwhm ./(2 .* (sqrt(2 .* log(2))));
% fwhm = 2 * sqrt(2*log(2)) * sigma

%% Start
clear all;close all;

N1 = lcfF2NErb(0.25); N2 = lcfF2NErb(8);
CF = 2; NOI = 15;

nRepCon = [1 10];
nStim = [4 8 16];
timeAcq = [3 7.5]; % TR in seconds
% loop modrate and sweep duration?
modRate = 1; % modulation rate of carrier - number of frequencies presented per second (Hz)
sweepDur = 2; % duration of sweep (seconds)
if modRate ==1
    stimType = 'discrete';
    nConPerSweep = 1;
else
    stimType = 'sweep';
    nConPerSweep = sweepDur * modRate; % number of conditions (frequencies presented) per sweep
    nRepCon = nRepCon * nConPerSweep;
end

for n = 1:length(nRepCon)
    nAcquisition(n,:) = nRepCon(n) .* nStim;
end
for i = 1:length(timeAcq)
    timeExp(:,:,i) = (nAcquisition .* timeAcq(i))./60; % experiemnet run time in minutes Row = number of repeats, Colum = number of Stimulus, Page = acquisition time
end
% for i = 1:length(sweepDur)
% nConPerSweep(i,:) = sweepDur(i) .* modRate; % number of conditions (frequencies presented) per sweep
% end
% make number of reps a function of total time aviable and conditions to
% test
% nStim = [4 64];
% acqTime = [3 7.5]; % TR in seconds
% expTime = 60; % experiement time in minutes
% nRep = round((expTime.*60) ./ acqTime);
% nRepCon = round(nRep / nStim .* nConPerSweep);
% create loop for sweep length
% create sweep vector - how to space?
% percent around or above CF or octave or ERB?
% spaced lin,log,erb,dlm
% devide repeats by number of conditions per sweep - times by images per
% scan
% stimFreq = nHz,CF,nStimID

C = [10 100];
C = C/100;
TW = CF * C;
npTW = length(C); % number of voxel Tuning widths to test
neCon = length(nStim); % number of stimulus sets to test
% fs = CF+(-3*CF:(6*CF)/(100-1):3*CF);
LowFreq = 0.1;
HighFreq = 10;
fs = LowFreq:0.01:HighFreq;

for n = 1:length(nRepCon)
    %% Predifine variables
    nRep = nRepCon(n);
    eCtrd = zeros(neCon,npTW); eSprd = zeros(neCon,npTW);
    dCtrd = zeros(neCon,npTW); dSprd = zeros(neCon,npTW);
    Ctrd = zeros(1,npTW); Sprd = zeros(1,npTW);
    fitCF = zeros(neCon,npTW); fitC = zeros(neCon,npTW); fitTW = zeros(neCon,npTW); fiterror = zeros(neCon,npTW);
    max_stimF = zeros(neCon,npTW);fwhm = zeros(neCon,npTW);
    
    for I = 1:neCon
        % create stimulus set frequencies
        stimF = lcfNErb2F(N1:((N2-N1)/(nStim(I)-1)):N2);
        if nConPerSweep == 1
            stimSweep = stimF;
        else
            for i = 1: length(stimF)
                stimSweep(:,i) = linspace(stimF(i),stimF(i)*2,nConPerSweep);
            end
            stimSweep = stimSweep/2;
        end
        for II = 1:npTW
            [Ct,Sp,f,w] = lcfW(CF,C(II));
            Ctrd(II) = Ct;
            Sprd(II) = Sp;
            jwd1 = zeros(1,nRep); jwd2 = zeros(1,nRep);
            jwd3 = zeros(1,nRep); jwd4 = zeros(1,nRep);
            %             r = zeros(1,nStim(I));
            r = zeros(size(stimSweep));
            for III = 1:nRep
                %                 [jwd1(III),jwd2(III),oneR] = lcfVoxResp(CF,C(II),NOI,stimF);
                [jwd1(III),jwd2(III),oneR] = lcfVoxResp(CF,C(II),NOI,stimSweep);
                jwd3(III) = abs(jwd1(III)-Ctrd(II))/Ctrd(II)*100;
                jwd4(III) = abs(jwd2(III)-Sprd(II))/Sprd(II)*100;
                r = r+oneR/nRep;
            end
            % need to seperate centroid function
            % bold response to responses
            eCtrd(I,II) = nanmean(jwd1); eSprd(I,II) = nanmean(jwd2);
            dCtrd(I,II) = nanmean(jwd3); dSprd(I,II) = nanmean(jwd4);
            
            xdata = stimSweep;
            if nConPerSweep == 1
                [max_value, max_index] = max(r);
                max_stimF(I,II) = stimSweep(max_index);
            else
                [max_r_value, max_r_index] = max(r);
                [max_value, max_index] = max(max_r_value);
                max_stimF(I,II) = stimSweep(max_r_index(max_index),max_index);
            end
            [fwhmIr, fwhmIc] = find(r>max_value/2);
            if numel(fwhmIr)== 1
                fwhm(I,II) = max_stimF(I,II) * 0.2;
            elseif numel(fwhmIr)== 0
                fwhm(I,II) = max_stimF(I,II) * 0.2;
            else
                % find highest and lowest value frequency
                % index all ones above threshold then find highest
                fwhm(I,II) = max(stimSweep(fwhmIr)) - min(stimSweep(fwhmIr));
                %                 fwhmIneg = [fwhmIr(1) fwhmIc(1)];
                %                 fwhmIpos = [fwhmIr(end) fwhmIc(end)];
                %                 fwhm(I,II) = abs(stimSweep(fwhmIpos(1),fwhmIpos(2)) - stimSweep(fwhmIneg(1),fwhmIneg(2)));
                %                 if fwhm(I,II) > (max_stimF(I,II)) *2
                %                     fwhm(I,II) = max_stimF(I,II);
                %                 end
            end
            fwhm(fwhm > max_stimF(I,II)) = (max_stimF(I,II));
            x0 = [max_stimF(I,II),fwhm(I,II)];
            ydata = r;
            %% think about this tomorrow
            if fwhm(I,II) < 0.01
                lb = [0.02, (max_stimF(I,II)*0.1)];
                x0 = [max_stimF(I,II),(max_stimF(I,II)*0.1)];
            else
                lb = [0.02, 0.01];
            end
            ub= [20, 100];
            fun = @lcfROEX;
            options = optimoptions('lsqcurvefit','TolFun',1e-10, 'MaxFunEvals', 400, 'MaxIter', 400); % default=1e-6
            [x, error] = lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options);
            
            fitCF(I,II) = x(1);
            fitC(I,II) = x(2);
            fitTW(I,II) = fitCF(I,II) .*fitC(I,II);
            fiterror(I,II)= error;
            fitW = lcfROEX (x,fs);
            efwhm(I,II) = 2 * sqrt(2*log(2)) * fitC(I,II);
            figure(n)
            %                 figure(length(nRepCon)+n)
            subplot(neCon,npTW,(neCon-I)*npTW+II), hold on
            for i = 1:length(stimF)
                plot(stimSweep(:,i),r(:,i),'ro--')
            end
            plot(f,w,'k-')
            plot(fs,fitW,'c--'),
            set(gca,'XScale','log')
            xlim ([LowFreq HighFreq])
            text(max(xlim),max(ylim),sprintf('Ctrd = %.2f (%.2f) kHz\nSprd = %.2f (%.2f) kHz\nCF = %.2f (%.2f) kHz\nTW = %.2f (%.2f) kHz\nx0 = %.2f, %.2f\nexpTimeCon = %.2f\nexpTimeSpar = %.2f',Ctrd(II),eCtrd(I,II),Sprd(II),eSprd(I,II),CF,fitCF(I,II),TW(II),fitTW(I,II),max_stimF(I,II),fwhm(I,II),timeExp(n,I,1),timeExp(n,I,2)),'HorizontalALignment','right','VerticalAlignment','top')
        end
    end
end

% figure(3), clf
% colormap('gray')
% subplot(1,2,1), contourf(C,NStim,dCtrd), set(gca,'XScale','log','YScale','log'), title('Ctrd'), xlabel('C'), ylabel('NStim')
% subplot(1,2,2), contourf(C,NStim,dSprd), set(gca,'XScale','log','YScale','log'), title('Sprd'), xlabel('C')
end
% Local functions
function nErb = lcfF2NErb(f)
nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);
end
function f = lcfNErb2F(nErb)
f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
end
function erb = lcfErb(f)
erb = 24.67*(4.37*f+1)/1000;
end
function [ECtrd,ESprd,r] = lcfVoxResp(CF,C,NOI,f)
% ERB = C*lcfErb(CF);
% P = 4*CF/ERB;
P = 4*CF/C;
g = abs(f-CF)/CF;
r = (1+P*g).*exp(-P*g)+NOI/100*randn(size(f));
ECtrd = sum(max(r,0).*f)/sum(max(r,0));
ESprd = sqrt(sum(max(r,0).*(f-ECtrd).^2)/sum(max(r,0)));
end
function [Ctrd,Sprd,f,w] = lcfW(CF,C)
% ERB = C*lcfErb(CF);
% P = 4*CF/ERB;
P = 4*CF/C;
f = CF+(-3*CF:(6*CF)/(100-1):3*CF); f = f(f>=0);
g = abs(f-CF)/CF;
w = (1+P*g).*exp(-P*g);

Ctrd = sum(w.*f)/sum(w);
Sprd = sqrt(sum(w.*(f-Ctrd).^2)/sum(w));
end
function [w] = lcfROEX(param,f)
CF = param(1);
C =param(2);
P = 4*CF/C;
% ERB = C*lcfErb(CF);
% P = 4*CF/ERB;
g = abs(f-CF)/CF;
w = (1+P*g).*exp(-P*g);
end
function hrf = makeHrf(TR)
% given the TR, return the HRF shape for t = 0 ... 30s
% using the equation given in the lecture (simple boynton version)
tau = 2; % time constant (s)
delta = 2; % time shift (s)
t = [0:TR:30]; % vector of time points (in steps of TR)
tshift = max(t-delta,0); % shifted, but not < 0
hrf = (tshift ./ tau) .^2 .* exp(-tshift ./tau )./(2*tau);
end