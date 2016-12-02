function runLSQSingle

% This function fits the Experiment-1 data; the data
% are adaptation in percent [(A-P)/A*100 = (1-P/A)*100, where A = unadapter response, P =
% adapted response as a function of the adapter-probe frequecy difference
% (df in cent) and adapter-probe SOA.
PF = 1; df = [0 200 600 1800]; af = PF*2.^(df/1200); 
soa = [125 250 500 1000]; Lev = 60;
y0 = funData(0,df); % These are the extrapolated data for SOA = 0 ms;

x = lsqnonlin(@(x)funSqESingle(PF,af,Lev,y0,x),[0.1 -10 1],[0.01 -25 0],[1 0 1]);
Ctip = 1; Ctail = x(1); W = x(2); Al = x(3);
save('LSQSingle.mat','Ctip','Ctail','W','Al')
% The parameters of the neuron frequency response functions are fitted to
% these data; lsqnonlin is part of the optimisation tool in Matlab; you may not need
% any fitting; 

ePProbe = funROEX(PF,Lev,Ctip,Ctail,W,false); PR = sum(ePProbe.excInt); % funROEX calculates the excitation pattern of a simulus (in this case, the probe);

figure(1), clf, hold on
plot(ePProbe.fc,ePProbe.excInt,'k'), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excInt') 

ad = zeros(1,length(af));
c = {'r' 'm' 'b' 'c'};
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,Ctip,Ctail,W,false); % excitation patter of the adapter;
    supp = 1-(ePAdapt.excInt/max(ePAdapt.excInt)).^Al; % this is the adaptation effect; it of proportional to a power-law function of eth adapter excitation.
    figure(1), plot(ePProbe.fc,supp.*ePProbe.excInt,'Color',c{I}), axis tight 
    ad(I) = (1-sum(supp.*ePProbe.excInt)/PR)*100; % this is the predicted adaptation [(1-P/A)*100].
end

figure(2), clf, hold on
plot(df,y0,'ks-','LineWidth',2)
plot(df,ad,'rs-','LineWidth',2)
text(max(xlim),max(ylim),sprintf('Ctip = 1, Ctail = %g, W = %g, Al = %g',Ctail,W,Al),'HorizontalAlignment','right','VerticalAlignment','top')
text(max(xlim),max(ylim)-0.05*diff(ylim),sprintf('Wtip = 1, Wtail = %g, Gw = %g',1/Ctail,10*log10(10^(W/10)/(1-10^(W/10)))),'HorizontalAlignment','right','VerticalAlignment','top')

ad = zeros(length(soa),length(af));
load('MeData.mat'); y = data(1:4,:);
for I = 1:length(af)
    ePAdapt = funROEX(af(I),Lev,Ctip,Ctail,W,false);
    for II = 1:length(soa)
        supp = 1-exp(-soa(II)/957.56)*(ePAdapt.excInt/max(ePAdapt.excInt)).^Al;
        ad(II,I) = (1-sum(supp.*ePProbe.excInt)/PR)*100;
    end
end
RMSD = sqrt(mean((ad(:)-y(:)).^2));

figure(3), clf, subplot(1,2,1), hold on
for I = 1:length(soa)
    plot(df,y(I,:),'ks-'), axis tight
    plot(df,ad(I,:),'ro--'), axis tight
end
set(gca,'YScale','log')
subplot(1,2,2), hold on
for I = 1:length(af)
    plot(soa,y(:,I),'ks-'), axis tight
    plot(soa,ad(:,I),'ro--'), axis tight
end
set(gca,'YScale','log')
text(max(xlim),max(ylim),sprintf('RMSD = %g%%',RMSD),'HorizontalAlignment','right','VerticalAlignment','top')

    
