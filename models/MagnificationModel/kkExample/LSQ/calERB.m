SF = 25/2; N = 2^15; DF = SF/N;
f = DF*(0:N-1); 
FC = 1; ERB1 = funErb(FC);
g = abs((f-FC)/FC);

Ptip = 4/ERB1;
Ptail = 0.0640389*4/ERB1;
W = 10^(-14.4381/10);
fws = (1-W)*(1+Ptip*g).*exp(-Ptip*g)+W*(1+Ptail*g).*exp(-Ptail*g);

Ptip = 1.8622*4/ERB1;
Ptail = 0.0941249*4/ERB1;
W = 10^(-23.3948/10);
fwm = (1-W)*(1+Ptip*g).*exp(-Ptip*g)+W*(1+Ptail*g).*exp(-Ptail*g);

ERBs = sum(fws*DF);
ERBm = sum(fwm*DF);
figure, hold on 
plot(f,fws), plot(f,fwm,'r')
text(max(xlim),max(ylim),sprintf('ERBs = %g, ERBm = %g, m/s = %g',ERBs,ERBm,ERBm/ERBs),'HorizontalAlignment','right','VerticalAlignment','top')
axis tight, set(gca,'XScale','log')


