function CriticalRatioPerERB_int = getCriticalRatioPerERB(f)

f = 10.^(linspace(log10(20),log10(20000),100));
CriticalRatio_int = zeros(1,length(f));

[CriticalRatio_int CriticalRatioMeasured f_measured_hz] = getCriticalRatio_HawkinsAndStevens1950(f);
ERB = getERB(f_measured_hz/1000);

CriticalRatioPerERB = CriticalRatioMeasured-10*log10(ERB);
CriticalRatio_Zwicker = 10*log10(10^(1/10)-1)*ones(size(f_measured_hz));

fit_extrapolate_100HzMinus = polyfit(f_measured_hz(1:4),CriticalRatioPerERB(1:4),1);
fit_extrapolate_9kHzPlus = polyfit(f_measured_hz(12:end),CriticalRatioPerERB(12:end),1);

CriticalRatioPerERB_int(f<100) = max(CriticalRatioPerERB(1),polyval(fit_extrapolate_100HzMinus,f(f<100)));
CriticalRatioPerERB_int(f>=100 & f<=9000) = interp1(f_measured_hz,CriticalRatioPerERB,f(f>=100 & f<=9000),'spline');
CriticalRatioPerERB_int(f>9000) = max(CriticalRatioPerERB(end),polyval(fit_extrapolate_9kHzPlus,f(f>9000)));

Nee_lev = -10*log10(ERB);
Nten_lev = -10*log10(ERB)-CriticalRatioPerERB;
Nee_CR_Zwicker_lev = -10*log10(ERB)-CriticalRatio_Zwicker;

figure; 
plot(f_measured_hz,CriticalRatioPerERB,'ko-');
hold on
plot(f,CriticalRatioPerERB_int,'r--')
legend('Critical Ratio (Measured)',...
    'Critical Ratio (Extrapolated)',...
    'Location','best')
set(gca,'XLim',[min(f) max(f)])
title('Critical Ratio Per ERB)')
xlabel('Frequency (Hz)') 
ylabel('dB')
end