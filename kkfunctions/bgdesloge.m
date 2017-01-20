function k = desloge

f = [0.125 0.25 0.5 1 1.5 2 3 4 6 8];
cr_1Hz = [17.75 16.3 17.25 18.5 19.25 20.5 22.5 25.1 26 27];
erb = 24.7*(4.37*f+1);
cr_ERB = cr_1Hz-10*log10(erb)

figure, clf, hold on
plot(f,cr_1Hz,'ko-')
plot(f,cr_ERB,'ro-')

% fft = logspace(0.02,20,10);
fft = 0.02:1:20;
cr_ERB_FR = interp1(f,cr_ERB,fft,'spline');

figure; hold on
semilogx(fft,cr_ERB_FR,'ko-')
semilogx(f,cr_ERB,'ro--')