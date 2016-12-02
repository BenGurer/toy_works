function data = LogVsLinModel
clf
% compF = [0.1]
% compF = [0.1 0.2 0.5 1 2 3 4 6 8]
a=  0.25;
b=  12;
n= 10;
% compF = logspace(log10(a),log10(b),n);
compF = linspace(a,b,n);
compLev= 1;
% compLev= ones(length(compF));
Ctip=1;
Ctail=1;
W=0;
PltFlag= 1;

for i = 1:length(compF)
data = funROEX(compF(i),compLev,Ctip,Ctail,W,PltFlag)
end

data.f = compF;
end

function excPat = funROEX(compF,compLev,Ctip,Ctail,W,PltFlag)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

excPat.nErb = funF2NErb(0.05):0.01:funF2NErb(20);
excPat.fc = funNErb2F(excPat.nErb);
excPat.erb = funErb(excPat.fc);
excPat.ptip = Ctip*4*excPat.fc./excPat.erb;
excPat.ptail = Ctail*4*excPat.fc./excPat.erb;
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
%     plot(excPat.fc,excPat.excLev,'k'), set(gca,'XScale','log'), xlabel('fc (kHz)'), ylabel('excLev (dB)')
    hold on
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
    end
    
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
