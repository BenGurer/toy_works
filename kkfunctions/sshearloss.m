function thr_SPL = sshearloss

ST = 1/25;
Nf = 2^14;
f = 1/(ST*2)*linspace(0,1,Nf);
Slope = 36.5; % dB/oct 

thr_HL = max(Slope*log2(f/3),0);
Thr_8 = thr_HL(find(f>=8,1,'first'));
thr_HL = min(thr_HL,Thr_8);

retSPL = xlsread('dBcalculator.xls','dB HL','A4:B14');
thr_SPL = thr_HL+interp1(retSPL(:,1)/1000,retSPL(:,2),f,'linear');

figure(2), clf
subplot(2,1,1), hold on
plot(retSPL(:,1)/1000,retSPL(:,2),'ro-')
title('Threshold of audibility ANSI S3.6 - 1996')
ylabel('dB SPL')
subplot(2,1,2), hold on
plot(f,thr_HL,'k--')
plot(f,thr_SPL,'b-')
legend('Thr in dB HL',...
    'Thr in dB SPL',...
    'Location','best')
set(gca,'XLim',[min(f) max(f)])
title('Steeply sloping hearig loss according to screening crtiterion')
xlabel('Frequency (kHz)') 
ylabel('dB (HL/SPL)')







