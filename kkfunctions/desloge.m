function k = desloge

f = [0.125 0.25 0.5 1 1.5 2 3 4 6 8];
cr_1Hz = [17.75 16.3 17.25 18.5 19.25 20.5 22.5 25.1 26 27];
erb = 24.7*(4.37*f+1);
cr_ERB = cr_1Hz-10*log10(erb);
cr_Zwicker = 10*log10(10^(1/10)-1)*ones(size(f));

figure, clf, hold on
plot(f,cr_1Hz,'ko-')
plot(f,cr_ERB,'ro-')
plot(f,cr_Zwicker,'b--')
legend('critical ratio (per Hz)',...
    'critical ratio (per ERB)',...
    'critical ratio (accroding to Zwicker)',...
    'Location','best')

