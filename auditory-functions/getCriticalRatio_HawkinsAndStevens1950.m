function CriticalRatio_int = getCriticalRatio_HawkinsAndStevens1950(f)
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
CriticalRatio = [19, 17.75, 17.5, 16.25, 16.5 17.25, 17.75, 18.5, 19.25, 20.5, 22.5, 25, 25.5, 26.75, 27, 28.5];
% CriticalRatioInt = interp1(f_measured_hz,CriticalRatio,f,'spline');

% f_extrapolate_100HzMinus = linspace(20,100,10);
% f_extrapolate_9kHzPlus = linspace(9000,20000,10); 
fit_extrapolate_100HzMinus = polyfit(f_measured_hz(1:4),CriticalRatio(1:4),1);
fit_extrapolate_9kHzPlus = polyfit(f_measured_hz(12:end),CriticalRatio(12:end),1);

CriticalRatio_int(f<100) = max(CriticalRatio(1),polyval(fit_extrapolate_100HzMinus,f(f<100)));
CriticalRatio_int(f>=100 & f<=9000) = interp1(f_measured_hz,CriticalRatio,f(f>=100 & f<=9000),'spline');
CriticalRatio_int(f>9000) = max(CriticalRatio(end),polyval(fit_extrapolate_9kHzPlus,f(f>9000)));

% [0.125 0.25 0.5 1 1.5 2 3 4 6 8];
% [17.75 16.3 17.25 18.5 19.25 20.5 22.5 25.1 26 27];

% semilogx(f_measured_hz,CriticalRatio); hold on
% semilogx(f,CriticalRatioInt,'-');

% plot(f,CriticalRatioInt,'-');hold on

figure; 
plot(f_measured_hz,CriticalRatio,'b');
hold on
% 
% plot(f_extrapolate_9kHzPlus,polyval(fit_extrapolate_9kHzPlus,f_extrapolate_9kHzPlus,'r--'));
% plot(f_extrapolate_100HzMinus,polyval(fit_extrapolate_100HzMinus,f_extrapolate_100HzMinus,'r--'));

plot(f,CriticalRatio_int,'k--')