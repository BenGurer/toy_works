function [CriticalRatio_int CriticalRatioMeasured f_measured_hz] = getCriticalRatio_HawkinsAndStevens1950(f)
%
%   usage: getCriticalRatio_HawkinsAndStevens1950
%      by: Ben Gurer
%    date: 19/01/2017
% purpose: output the critical ratio values as determined by Hawkins And Stevens (1950)
% 
% discription:
% Measured values taken from Hawkins And Stevens (1950) The masking of pure tones and of speech by white noise (fig. 6)
% "Ratio between the monaural masked threshold of a pure tone and the level per cycle of the masking noise measured at the frequency of the pure tone".

CriticalRatio_int = zeros(1,length(f));

f_measured_hz = [100, 125, 175, 250, 350, 500, 700, 1000, 1400, 2000, 2800, 4000, 5600, 7000, 8000, 9000];
CriticalRatioMeasured = [19, 17.75, 17.5, 16.25, 16.5 17.25, 17.75, 18.5, 19.25, 20.5, 22.5, 25, 25.5, 26.75, 27, 28.5];
% CriticalRatioInt = interp1(f_measured_hz,CriticalRatio,f,'spline');

% f_extrapolate_100HzMinus = linspace(20,100,10);
% f_extrapolate_9kHzPlus = linspace(9000,20000,10); 
fit_extrapolate_100HzMinus = polyfit(f_measured_hz(1:4),CriticalRatioMeasured(1:4),1);
fit_extrapolate_9kHzPlus = polyfit(f_measured_hz(12:end),CriticalRatioMeasured(12:end),1);

CriticalRatio_int(f<100) = max(CriticalRatioMeasured(1),polyval(fit_extrapolate_100HzMinus,f(f<100)));
CriticalRatio_int(f>=100 & f<=9000) = interp1(f_measured_hz,CriticalRatioMeasured,f(f>=100 & f<=9000),'spline');
CriticalRatio_int(f>9000) = max(CriticalRatioMeasured(end),polyval(fit_extrapolate_9kHzPlus,f(f>9000)));


% figure; 
% plot(f_measured_hz,CriticalRatioMeasured,'ko-');
% hold on
% plot(f,CriticalRatio_int,'r--')
% legend('Critical Ratio (Measured)',...
%     'Critical Ratio (Extrapolated)',...
%     'Location','best')
% set(gca,'XLim',[min(f) max(f)])
% title('Critical Ratio - Hawkins And Stevens (1950)')
% xlabel('Frequency (kHz)') 
% ylabel('dB')
end