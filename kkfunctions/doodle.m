function doodle

FS = 25;
N = 2*2^nextpow2(FS/2*1000);

% a time signal with length N and rms = 1
noi = randn(1,N);
noi = noi/(norm(noi)/sqrt(N));
sqrt(mean(noi.^2))

% the rms of the fft is then sqrt(N), not 1
frq = FS/2*linspace(0,1,N/2);
NOI = fft(noi); 
sqrt(mean(abs(NOI).^2)), sqrt(N)

% but converted back, the rms is still 1
nnoi = real(ifft(NOI));
sqrt(mean(nnoi.^2))

% if the signalis filtered on the way, the resulting rms will be equal to
% the  rms of the filter 
bp = zeros(size(frq)); 
bp(and(frq>=funNErb2F(funF2NErb(1)-0.5),frq<=funNErb2F(funF2NErb(1)+0.5))) = 1;
% bp = bp/sqrt(mean(bp.^2));
sqrt(mean(abs(bp).^2))
nnoi = real(ifft([bp fliplr(bp)].*NOI));
sqrt(mean(nnoi.^2))

% thus if we scale the ee filter, such as its ERB-rms is 1 ...
ee = 1./sqrt(funErb(frq));
% ee = ee/sqrt(mean(ee.^2));
% sqrt(mean(abs(ee).^2))
% nnoi = real(ifft([ee fliplr(ee)].*NOI));
% sqrt(mean(nnoi.^2))

% ... then we get an appropriately scaled time signal, obviating eth need
% to scale with LEE!
ee = ee/sqrt(mean((bp.*ee).^2));
sqrt(mean(abs(ee).^2)), 20*log10(sqrt(mean(abs(ee).^2)))
nnoi = real(ifft([ee fliplr(ee)].*NOI));
sqrt(mean(nnoi.^2)), 20*log10(sqrt(mean(nnoi.^2)))






