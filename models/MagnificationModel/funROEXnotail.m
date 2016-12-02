function excPat = funROEX(compF,compLev,Ctip,PltFlag)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);
n=4;
excPat.nErb = funF2NErb(0.05):0.01:funF2NErb(15);
excPat.fc = funNErb2F(excPat.nErb);
excPat.erb = funErb(excPat.fc);
excPat.ptip = Ctip*4*excPat.fc./excPat.erb;
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
        fw=(1+(excPat.ptip(I)*g)).*exp(-excPat.ptip(I)*g); %Two parameter roex
%         fw=exp(-excPat.ptip(I)*g); %Exponential
%         fw = k*abs(((1 + 1i*(excPat.ptip(I)-1)/b).^-n + (1 + 1i*(excPat.ptip(I)+1)/b).^-n));
        excPat.excInt(I) = excPat.excInt(I)+fw*compInt(II);
    end
end
excPat.excLev = 10*log10(max(excPat.excInt,1));

if PltFlag==true
    figure(1), clf
    plot(excPat.fc,excPat.excLev,'k'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)')
    for I = 1:length(compF)
        line(repmat(compF(I),1,2),[min(ylim) compLev(I)],'Color','r')
    end
end


