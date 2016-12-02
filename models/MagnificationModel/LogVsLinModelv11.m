function data = LogVsLinModel
figure(1), clf

StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10;
% compF = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps)
compF = linspace(StimLowFreq,StimHighFreq,StimSteps)
% compF = logspace(log10(0.5),log10(2),100);
% compF = 1
% compLev= 1;
compLev= ones(length(compF));
C=0.06;
Ctail=0;
W=0;
PltFlag= 0;
plotC = {'k','b','r','g','y'}; % Cell array of colros.
plotC = repmat(plotC,1,round(length(compF))./5);
for i = 1:length(compF)
    excPat  = funROEX(compF(i),compLev,C,PltFlag);
    
    figure(1)
    plot(excPat.fc,excPat.excLev,'Color',plotC{i}), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)')
    hold on
%     for I = 1:length(compF)
        line(repmat(compF(i),1,2),[min(ylim) compLev(i)],'Color',plotC{i})
%         hold on
%     end
    
end
% excPat = funROEX(compF,compLev,C,PltFlag);

end

function excPat = funROEX(compF,compLev,C,PltFlag)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

HighFreq = 20;
LowFreq = 0.05;
NeuronSpacing = ((HighFreq - LowFreq)./100).*1;
% C = ((HighFreq - LowFreq)./100).*1
% NeuronSpacing = 400
% NeuronSpacing = 0.01
excPat.nErb = funF2NErb(LowFreq):NeuronSpacing:funF2NErb(HighFreq);
% excPat.nErb = linspace(funF2NErb(LowFreq),funF2NErb(HighFreq),NeuronSpacing);
% excPat.nErb = logspace(log10(funF2NErb(LowFreq)),log10(funF2NErb(HighFreq)),NeuronSpacing);
excPat.fc = funNErb2F(excPat.nErb);
excPat.erb = funErb(excPat.fc);
excPat.p = C*4*excPat.fc./excPat.erb;
% excPat.ptail = Ctail*4*excPat.fc./excPat.erb;
% if isinf(W)
%     W = 0;
% else
%     W = 10^(W/10);
% end

excPat.excInt = zeros(1,length(excPat.fc));
for I = 1:length(excPat.fc)
    for II = 1:length(compF)
        g = abs((compF(II)-excPat.fc(I))/excPat.fc(I));
        % fw = (1-W)*(1+excPat.ptip(I)*g).*exp(-excPat.ptip(I)*g)+W*(1+excPat.ptail(I)*g).*exp(-excPat.ptail(I)*g);
        fw=(1+(excPat.p(I)*g)).*exp(-excPat.p(I)*g); %Two parameter roex
        %         fw=exp(-excPat.ptip(I)*g); %Exponential
        %         fw = k*abs(((1 + 1i*(excPat.ptip(I)-1)/b).^-n + (1 + 1i*(excPat.ptip(I)+1)/b).^-n));
        excPat.excInt(I) = excPat.excInt(I)+fw*compInt(II);
    end
end
excPat.excLev = 10*log10(max(excPat.excInt,1));
if PltFlag==true
    figure(1)
    plot(excPat.fc,excPat.excLev,'g'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)')
    hold on
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
        hold on
    end
end

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
