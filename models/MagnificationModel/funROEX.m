function excPat = funROEX(compF,compLev,Ctip,Ctail,W,PltFlag)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are the widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

excPat.nErb = funF2NErb(0.05):0.01:funF2NErb(15);  
excPat.fc = funNErb2F(excPat.nErb); % ERB-scaled filter frequency axis;
excPat.erb = funErb(excPat.fc); % filter widths according to G&M 1990;
excPat.ptip = Ctip*4*excPat.fc./excPat.erb; % resulting p values for tip and tail; 
excPat.ptail = Ctail*4*excPat.fc./excPat.erb;
if isinf(W)
    W = 0;
else
    W = 10^(W/10);
end

excPat.excInt = zeros(1,length(excPat.fc));
for I = 1:length(excPat.fc) 
    for II = 1:length(compF)
        G = abs((compF(II)-excPat.fc(I))/excPat.fc(I)); % relative frequency for fc(I);
        FW = (1-W)*(1+excPat.ptip(I)*G).*exp(-excPat.ptip(I)*G)+W*(1+excPat.ptail(I)*G).*exp(-excPat.ptail(I)*G); % filter weights for fc(I) at compF(II); 
        excPat.excInt(I) = excPat.excInt(I)+FW*compInt(II); % filter response for fc(I);
    end
end
excPat.excLev = 10*log10(max(excPat.excInt,1)); % convert filter responses to decibels; 

if PltFlag==true
    figure(1), clf
    plot(excPat.fc,excPat.excLev,'k'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)') 
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
    end
end


