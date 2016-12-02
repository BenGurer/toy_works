function data = LogVsLinModel
clf

HighFreq = 20;
LowFreq = 0.05;
NeuronSpacing = ((HighFreq - LowFreq)./100).*10;

StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10
% compF = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps);
% compF = linspace(StimLowFreq,StimHighFreq,StimSteps)
compF = 1
compLev= 1;
Ctip=0.1;
Ctail=0.1;
W=0;
PltFlag= 1;

for i = 1:length(compF)
data = funROEX(compF(i),compLev,Ctip,Ctail,W,PltFlag)
end

data.f = compF;
data.Nspread = NeuronSpacing; 
end

function excPat = funROEX(compF,compLev,Ctip,Ctail,W,PltFlag)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

HighFreq = 20;
LowFreq = 0.05;
% NeuronSpacing = ((HighFreq - LowFreq)./100).*10;
NeuronSpacing = 0.001
excPat.nErb = funF2NErb(LowFreq):NeuronSpacing:funF2NErb(HighFreq);
excPat.fc = funNErb2F(excPat.nErb);
excPat.erb = funErb(excPat.fc);
excPat.ptip = Ctip*4*excPat.fc./excPat.erb;
excPat.ptail = Ctail*4*excPat.fc./excPat.erb;
sigma = NeuronSpacing./100;

if isinf(W)
    W = 0;
else
    W = 10^(W/10);
end

excPat.excInt = zeros(1,length(excPat.fc));
for I = 1:length(excPat.fc) 
    for II = 1:length(compF)
        g = abs((compF(II)-excPat.fc(I))/excPat.fc(I));
        fw = (1-W)*(1+excPat.ptip(I)*g).*exp(-excPat.ptip(I)*g)+W*(1+excPat.ptail(I)*g).*exp(-excPat.ptail(I)*g);
        excPat.excInt(I) = excPat.excInt(I)+fw*compInt(II);
    end
end
excPat.excLev = 10*log10(max(excPat.excInt,1));

if PltFlag==true
%     figure(1), clf

    figure(1)
    plot(excPat.fc,excPat.excLev,'k'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)') 
    hold on
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
    end
    figure(2)
    plot(excPat.fc,fw), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)') 
end
end
function erb = funErbGM90(f)
% Erb (in kHz) according to Glasberg and Moore 1990 as a function of frequency, f (in kHz); 
% f can be a vector.

erb = 24.67*(4.37*f+1)/1000;
end
function f = funNErb2FGM90(nErb)
% Converts Erb number (nErb) to frequency, f (in kHz); nErb can be a vector.

f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
end
function nErb = funF2NErbGM90(f)
% Converts frequency, f (in kHz), to Erb number (nErb) according to Glasberg and Moore 1990; 
% f can be a vector; nErb = integral from f = 0 to f of 1/erb(f).

nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);
end
