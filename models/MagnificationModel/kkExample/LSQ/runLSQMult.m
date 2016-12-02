function runLSQMult

PF = 1; df = [0 200 600 1800]; af = PF*2.^(df/1200); Lev = 60;

load('MeData.mat'); y = data(5:6,:); y0 = exp(mean(log(y)));
load('LSQSingle.mat'); M = 1.044321/1000;
ERBs = funCalERB(Ctip,Ctail,W); SnglCtail = Ctail;

x = lsqnonlin(@(x)funSqEMult(PF,af,Lev,Al,M,y0,x),[1 0.1 -15],[0.1 0.01 -50],[10 1 0]);
Ctip = x(1); Ctail = x(2); W = x(3); 
save('LSQMult.mat','Ctip','Ctail','W','Al')
ERBm = funCalERB(Ctip,Ctail,W); MultCtail = Ctail;

ePProbe = funROEX(PF,Lev,x(1),x(2),x(3),false);
PR = sum(ePProbe.excInt);

figure(4), clf, hold on
plot(ePProbe.fc,ePProbe.excInt,'k'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excInt') 

ad2 = zeros(1,length(af));
ad3 = zeros(1,length(af));
ad = zeros(1,length(af));
c = {'r' 'm' 'b' 'c'};
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,Ctip,Ctail,W,false);
    supp2 = (1-exp(-M*1000)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al);
    supp3 = (1-exp(-M*1500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*1000)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al).*...
        (1-exp(-M*500)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al);
    figure(4), plot(ePProbe.fc,mean([supp2;supp3],1).*ePProbe.excInt,'Color',c{I}), axis tight 
    ad2(I) = (1-sum(supp2.*ePProbe.excInt)/PR)*100;
    ad3(I) = (1-sum(supp3.*ePProbe.excInt)/PR)*100;
    ad(I) = exp(mean(log([ad2(I) ad3(I)])));
end

figure(5), clf, hold on
plot(df,y0,'ks-','LineWidth',2)
plot(df,ad,'rs-','LineWidth',2)
text(max(xlim),max(ylim),sprintf('Ctip = %g, Ctail = %g, W = %g, Al = %g',Ctip,Ctail,W,Al),'HorizontalAlignment','right','VerticalAlignment','top')
text(max(xlim),max(ylim)-0.05*diff(ylim),sprintf('Wtip = %g, Wtail = %g, Gw = %g',1/Ctip,1/Ctail,10*log10(10^(W/10)/(1-10^(W/10)))),'HorizontalAlignment','right','VerticalAlignment','top')
text(max(xlim),max(ylim)-0.1*diff(ylim),sprintf('Ctail m/s = %g, ERB s/m = %g',MultCtail/SnglCtail,ERBs/ERBm),'HorizontalAlignment','right','VerticalAlignment','top')

RMSD = sqrt(mean(([ad2 ad3]-[y(1,:) y(2,:)]).^2));

figure(6), clf, hold on
plot(df,y(1,:),'ks-'), plot(df,y(2,:),'ko--')
plot(df,ad2,'rs-'), plot(df,ad3,'ro--')
set(gca,'YScale','log')
text(max(xlim),max(ylim),sprintf('RMSD = %g%%',RMSD),'HorizontalAlignment','right','VerticalAlignment','top')


