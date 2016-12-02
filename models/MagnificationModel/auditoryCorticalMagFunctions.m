f=(200:10:8000);

dlf = @(f)10.^(0.026*sqrt(f)-0.533);

erb = @(f)24.7*(4.37*f/1000+1);


figure;plot(f,dlf(f));

figure;plot(f,dlf(f)./f);

figure;plot(f,erb(f));

figure;plot(f,erb(f)./f);


figure;plot(sqrt(f),log10(dlf(f)));
set(gca,'xlim',sqrt([0 10000]),'xTick',sqrt([0 0.2 0.4 0.6 0.8 1 2 4 8]*1000));
set(gca,'xTickLabel',num2str((get(gca,'xTick').^2/1000)'));
set(gca,'ylim',log10([0.1 5000]),'yTick',log10([0.5 1 5 10 50 100 500 1000]));
set(gca,'yTickLabel',num2str((10.^(get(gca,'yTick')))'));