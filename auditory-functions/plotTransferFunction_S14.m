function tf = plotTransferFunction_S14
filename{1}= ('N:/matlab/tdtMRI/transferFunctions/S14_396insertsLeftFFT.csv');
filename{2}= ('N:/matlab/tdtMRI/transferFunctions/S14_396insertsRightFFT.csv');
tf = struct([]);
figure
for i = 1:length(filename)
[path,file,extension] = fileparts(filename{i});
switch(extension)
  case '.csv'  % transfer functions measured using IHR Brüel&Kjær microphone and 2cc coupler
    ffts = csvread(filename{i},29,1);
    ffts = ffts(:,1:end-3); %remove last 3 columns corresponding to two different weighted averages and an empty column

    tf(i).frequencies = (0:size(ffts,2)-1)*(20000/(size(ffts,2)-1)); % Frequency step size based on number of sample points
    tf(i).frequencies = tf(i).frequencies/1000; %convert to kHz
    tf(i).fft = mean(ffts(6:end,:));
    tf(i).fft = conv(tf(i).fft,ones(1,20)/20,'same');    
    
  case '.bin' %Sensimetrics S14 impulse response provided by manufacturer
    [h,Fs] = load_filter(filename{i});
    impulseResponse=h;
    zeroLocation = 25; %set arbitrary 0 time sample just before impulse reponse
    Fs = Fs/1000;%convert to ms
    time = ((1:length(h))' - zeroLocation)*1/Fs; %convert to ms
    impulse = zeros(size(time));
    impulse(zeroLocation)= 0.25; %pulse with arbitrary voltage

    % compute FFT
    nFFT = 2^nextpow2(size(time,1));
    transferFreqResolution = Fs/nFFT;
    tf(i).frequencies = (0:nFFT/2-1)*transferFreqResolution;

    impulseResponseFft = 20*log10(abs(fft(impulseResponse,nFFT)));
    impulseFft = 20*log10(abs(fft(impulse,nFFT)));
    impulseResponseFft = impulseResponseFft(1:end/2);
    impulseFft = impulseFft(1:end/2);
    tf(i).fft = impulseFft - impulseResponseFft;  %not sure why I need to take the negative of the impulse response
end
%centre on 1kHz
[~,f1kHz] = min(abs(tf(i).frequencies-1));
tf(i).fft = tf(i).fft - tf(i).fft(f1kHz);

plot(tf(i).frequencies,-(tf(i).fft))
hold on
end

% %cap at minAttenuation dB
% minAttenuation = -20;
% tf.fft(tf.fft<minAttenuation) = minAttenuation;

set(gca,'XLim',[min(tf(i).frequencies) max(tf(i).frequencies)])
title('S14 Transfer Function (normalised at 1kHz)')
xlabel('Frequency (kHz)') 
ylabel('dB SPL')
legend('Left',...
    'Right',...
    'Location','best')


