function runLSQSharp

PF = 1; df = [0 200 600 1800]; af = PF*2.^(df/1200); 
soa = [125 250 500 1000]; Lev = 60;
y0 = funData(0,df);

x = lsqnonlin(@(x)funSqESharp(PF,af,Lev,y0,x),[0.1 -10 1 5],[0.01 -25 0 1],[1 0 1 15]);
Ctip = 1; Ctail = x(1); W = x(2); Al = x(3); S = x(4);

ePProbe = funROEX(PF,Lev,1,Ctail,W,false);
PR = sum(ePProbe.excInt);

figure(6), clf, hold on
plot(ePProbe.fc,ePProbe.excInt,'k'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excInt') 

ad = zeros(1,length(af));
c = {'r' 'm' 'b' 'c'};
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,Ctip,Ctail,W,false);
    ePProbe = funROEXSharp(PF,Lev,Ctip,Ctail,W,(ePAdapt.excInt/max(ePAdapt.excInt)).^Al,S,false);
    figure(6), plot(ePProbe.fc,ePProbe.excInt,'Color',c{I}), axis tight 
    ad(I) = (1-sum(ePProbe.excInt)/PR)*100;
end

figure(7), clf, hold on
plot(df,y0,'ks-','LineWidth',2)
plot(df,ad,'rs-','LineWidth',2)
text(max(xlim),max(ylim),sprintf('Ctip = 1, Ctail = %g, W = %g, Al = %g, S = %g',Ctail,W,Al,S),'HorizontalAlignment','right','VerticalAlignment','top')
text(max(xlim),max(ylim)-0.05*diff(ylim),sprintf('Wtip = 1, Wtail = %g, Gw = %g',1/Ctail,10*log10(10^(W/10)/(1-10^(W/10)))),'HorizontalAlignment','right','VerticalAlignment','top')

s = zeros(1,length(soa));
ad = zeros(length(soa),length(af));
load('MeData.mat'); y = data(1:4,:);
for I = 1:length(soa)
    s(I) = lsqnonlin(@(x)funSqESOA(PF,af,Lev,Ctail,W,Al,y(I,:),x),5,1,15);
    for II = 1:length(af)
        ePAdapt = funROEX(af(II),Lev,Ctip,Ctail,W,false);
        ePProbe = funROEXSharp(PF,Lev,Ctip,Ctail,W,(ePAdapt.excInt/max(ePAdapt.excInt)).^Al,s(I),false);
        ad(I,II) = (1-sum(ePProbe.excInt)/PR)*100;
    end
end
text(max(xlim),max(ylim)-0.1*diff(ylim),sprintf('s = %g, %g, %g, %g',s(1),s(2),s(3),s(4)),'HorizontalAlignment','right','VerticalAlignment','top')
RMSD = sqrt(mean((ad(:)-y(:)).^2));

figure(8), clf, subplot(1,2,1), hold on
for I = 1:length(soa)
    plot(df,y(I,:),'ks-');
    plot(df,ad(I,:),'ro--'); hold on
end
set(gca,'YScale','log')

subplot(1,2,2), hold on
for I = 1:length(af)
    plot(soa,y(:,I),'ks-'), axis tight
    plot(soa,ad(:,I),'ro--'), axis tight
end
set(gca,'YScale','log')
text(max(xlim),max(ylim),sprintf('RMSD = %g%%',RMSD),'HorizontalAlignment','right','VerticalAlignment','top')

figure(9), clf, hold on
bar(soa,s,'FaceColor',0.5*ones(1,3)), plot(soa,exp(-soa/957.56)*max(s),'ks-')
set(gca,'YScale','log')

save('LSQSharp.mat','Ctip','Ctail','W','Al','s')




