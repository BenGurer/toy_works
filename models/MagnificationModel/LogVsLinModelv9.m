function data = LogVsLinModel
    figure(1), clf

HighFreq = 20;
LowFreq = 0.05;
NeuronSpacing = ((HighFreq - LowFreq)./100).*10;

StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10;
% compF = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps);
% compF = linspace(StimLowFreq,StimHighFreq,StimSteps)
compF = logspace(log10(0.5),log10(2),100);
% compF = 1
% compLev= 1;
compLev= ones(length(compF));
C=1;
Ctail=0;
W=0;
PltFlag= 1;

% for i = 1:length(compF)
% data = funROEX(compF(i),compLev,Ctip,Ctail,W,PltFlag)
% end
data = funROEX(compF,compLev,C,PltFlag)
data.f = compF;
data.Nspread = NeuronSpacing; 
end

function excPat = funROEX(compF,compLev,C,PltFlag)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

HighFreq = 20;
LowFreq = 0.05;
% NeuronSpacing = ((HighFreq - LowFreq)./100).*10;
% NeuronSpacing = 400
NeuronSpacing = 0.1
excPat.nErb = funF2NErb(LowFreq):NeuronSpacing:funF2NErb(HighFreq);
% excPat.nErb = linspace(funF2NErb(LowFreq),funF2NErb(HighFreq),NeuronSpacing);
% excPat.nErb = logspace(log10(funF2NErb(LowFreq)),log10(funF2NErb(HighFreq)),NeuronSpacing);
excPat.fc = funNErb2F(excPat.nErb);
excPat.erb = funErb(excPat.fc);
% excPat.ptip = Ctip*4*excPat.fc./excPat.erb;
% excPat.ptail = Ctail*4*excPat.fc./excPat.erb;
% p = Ctip*4*excPat.fc./excPat.erb;
p = C*4*excPat.fc./excPat.erb;
% p = Ctail*4*excPat.fc./excPat.erb;
cF = 139;
excPat.cF = excPat.fc(cF);
% aCf = 1;

% if isinf(W)
%     W = 0;
% else
%     W = 10^(W/10);
% end

excPat.excInt = zeros(1,length(excPat.fc));
excPat.FilRes = zeros(1,length(compF));
% excPat.FilRes = zeros(10);
% for I = 1:length(excPat.fc) 
    for II = 1:length(compF)
%         g = abs((compF(II)-excPat.fc(I))/excPat.fc(I));
%         g = abs((compF(II)-aCf)/aCf);
        g = abs((compF(II)-excPat.cF)/excPat.cF);
        fw = (1+(p(cF))*g).*exp((-p(cF))*g);
%         fw = (1+(p)*g).*exp((-p))*g;
        excPat.FilRes(II) = excPat.FilRes(II)+fw*compInt(II);
%         excPat.excInt(I) = excPat.excInt(I)+fw;
%         excPat.excInt(I) = excPat.excInt(I)+fw*compInt(II);
%         excPat.excInt(I) = fw*compInt(II);
    end
% end
excPat.excLev = 10*log10(max(excPat.FilRes,1));


if PltFlag==true
    figure(1)
    plot(compF,excPat.excLev ,'g'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)') 
    hold on
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
        hold on
    end
 end  

% if PltFlag==true
%     figure(1)
%     plot(excPat.fc,excPat.excLev,'g'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)') 
%     hold on
%     for I = 1:length(compF)
%         line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
%         hold on
%     end
%  end  

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
