function transformingvectorspace
xmin = 1;
xmax = 100;
steps = 100;

xLin = linspace(xmin,xmax,steps);

xLog = 10.^(xLin);

dataLog = logspace(log10(xmin),log10(xmax),steps);

plot (xLog,dataLog)

f0=1000;
y=(2.^([0:100:3000]./1200)).*f0;

end