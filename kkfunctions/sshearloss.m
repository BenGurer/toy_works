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

[f_hz, ThresholdOfHearing_FreeField_db_SPL, ThresholdOfHearing_DiffuseField_db_SPL] = getThresholdOfHearing;
[f_measured_hz, CriticalRatio] = getCriticalRatio_HawkinsAndStevens1950;
% W = 2*pi*f;
erb = 24.7*(4.37*f+1);
erb = 24.7 * (4.37 * f / 1000 + 1);
k_bells = 10*log10(2*pi*erb);


figure, clf
subplot(2,1,1), hold on
semilogx(retSPL(:,1),retSPL(:,2),'ro-')
semilogx(f_hz,ThresholdOfHearing_DiffuseField_db_SPL)
semilogx(f_hz,ThresholdOfHearing_FreeField_db_SPL)
semilogx(f_measured_hz,CriticalRatio)
plot(f,k_bells,'k--')
% plot(f,W,'b--')
plot(f,erb,'g--')
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







