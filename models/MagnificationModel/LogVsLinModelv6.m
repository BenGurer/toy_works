function data = LogVsLinModel
close all; clc

HighFreq = 20;
LowFreq = 0.05;
NeuronSpacing = ((HighFreq - LowFreq)./100).*10;

StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10
% compF = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps);
% compF = linspace(StimLowFreq,StimHighFreq,StimSteps)
compF = [linspace(0.5,2,10)]
% compF = 1
% compLev= 1;
compLev= ones(length(compF));
Ctip=0.06;
Ctail=0;
W=0;
PltFlag= 1;

% for i = 1:length(compF)
% data = funROEX(compF(i),compLev,Ctip,Ctail,W,PltFlag)
% end
data = funROEX(compF,compLev,Ctip,Ctail,W,PltFlag)
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
NeuronSpacing = 0.1
excPat.nErb = funF2NErb(LowFreq):NeuronSpacing:funF2NErb(HighFreq);
% excPat.nErb = funF2NErb(1);
excPat.fc = funNErb2F(excPat.nErb);
excPat.erb = funErb(excPat.fc);
% excPat.ptip = Ctip*4*excPat.fc./excPat.erb;
% excPat.ptail = Ctail*4*excPat.fc./excPat.erb;
p = Ctip*4*excPat.fc./excPat.erb;
% p = Ctail*4*excPat.fc./excPat.erb;
cF = 139;
% excPat.ptipkH = zeros(1,length(excPat.ptip));
% excPat.ptipkH(cF) = excPat.ptip(cF);
% excPat.ptailkH = zeros(1,length(excPat.ptail));
% excPat.ptailkH(cF) = excPat.ptail(cF);
% sigma = NeuronSpacing./100;

if isinf(W)
    W = 0;
else
    W = 10^(W/10);
end

excPat.excInt = zeros(1,length(excPat.fc));
for I = 1:length(excPat.fc) 
    for II = 1:length(compF)
        g = abs((compF(II)-excPat.fc(cF))/excPat.fc(cF));
        fw = (1+p(cF)*g).*exp(-p(cF)*g);
        excPat.excInt(I) = excPat.excInt(I)+fw*compInt(II);
    end
end
excPat.excLev = 10*log10(max(excPat.excInt,1));

if PltFlag==true
%     figure(1), clf

    figure(1)
    plot(excPat.fc,excPat.excLev,'g'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)') 
    hold on
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
        hold on
    end
 end  
%     cF = 139;
%     
% excPat.FilRes = zeros(1,length(excPat.fc));
% for I = 1:length(excPat.fc) 
%     for II = 1:length(compF)
%         g = abs((compF(II)-excPat.fc(I))/excPat.fc(I));
%         fw = (1-W)*(1+excPat.ptipkH(II)*g).*exp(-excPat.ptipkH(II)*g)+W*(1+excPat.ptailkH(II)*g).*exp(-excPat.ptailkH(II)*g);
%         excPat.FilRes(I) = excPat.FilRes(I)+fw*compInt(II);
%     end
% end
% excPat.excLev = 10*log10(max(excPat.FilRes,1));
%     figure(1)
%     plot(excPat.fc,excPat.excLev, 'g')
% figure(2)
% plot(excPat.fc,excPat.ptip, 'r'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)')
end
function erb = funErb(f)
% Erb (in kHz) according to Glasberg and Moore 1990 as a function of frequency, f (in kHz); 
% f can be a vector.

erb = 24.67*(4.37*f+1)/1000;
end
function f = funNErb2F(nErb)
% Converts Erb number (nErb) to frequency, f (in kHz); nErb can be a vector.

f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
end
function nErb = funF2NErb(f)
% Converts frequency, f (in kHz), to Erb number (nErb) according to Glasberg and Moore 1990; 
% f can be a vector; nErb = integral from f = 0 to f of 1/erb(f).

nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);
end
