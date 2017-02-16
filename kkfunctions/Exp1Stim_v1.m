function Exp1Stim

Fl = 0.25; Fh = 8; % Fl/h = low/high cutoff freuqnecy;

ST = 1/50;
D = 0.4356;
K = 0.0920*1000^(1-D); 

Dur = 10;
N = round(Dur/ST); tim = (0:(N-1))*ST; 

TWin = 1.25;
twin = [cos(pi/2*(tim(find(tim<=TWin))-TWin)/TWin).^2 ...
    ones(1,length(find(and(tim>TWin,tim<Dur-TWin)))) ...
    cos(pi/2*(tim(find(tim>=Dur-TWin))-Dur+TWin)/TWin).^2]; 

DF = 1/100; % 10-Hz freuqency sampling;
Ml = round(Fl/DF); Mh = round(Fh/DF);
fl = lcfNErb2F(lcfF2NErb(Fl)+[-0.5 0.5]); fh = lcfNErb2F(lcfF2NErb(Fh)+[-0.5 0.5]);
ml = [ceil(fl(1)/DF) floor(fl(2)/DF)]; mh = [ceil(fh(1)/DF) floor(fh(2)/DF)]; 
fwin = [cos(pi/2*((ml(1):ml(2))*DF-fl(2))/diff(fl)) ones(1,mh(1)-ml(2)-1) cos(pi/2*((mh(1):mh(2))*DF-fh(1))/diff(fh))]; % spectral window;

aw = ones(1,mh(2)-ml(1)+1).*fwin; aw = aw*sqrt(2)/sum(aw);
ap = 1./sqrt(lcfERB((ml(1):mh(2))*DF)).*fwin; ap = ap*sqrt(sum(aw.^2)/sum(ap.^2));
phi = lcfPhase((ml(1):mh(2))*DF,K,D);

Nt = 100*10/ST;
Nfft = 2^nextpow2(Nt);
frq = 1/(ST*2)*linspace(0,1,Nfft/2);

clkw = zeros(size(tim));
clkp = zeros(size(tim));
chpw = zeros(size(tim));
chpp = zeros(size(tim));
for I = ml(1):mh(2)
    clkw = clkw+aw(I-ml(1)+1)*cos(2*pi*I*DF*(tim-Dur/2)); 
    clkp = clkp+ap(I-ml(1)+1)*cos(2*pi*I*DF*(tim-Dur/2)); 
    chpw = chpw+aw(I-ml(1)+1)*cos(2*pi*I*DF*(tim-Dur)-phi(I-ml(1)+1)); 
    chpp = chpp+ap(I-ml(1)+1)*cos(2*pi*I*DF*(tim-Dur)-phi(I-ml(1)+1)); 
end
clkw = clkw.*twin;
clkp = clkp.*twin;
chpw = chpw.*twin;
chpp = chpp.*twin;

spclkw = fft(clkw,Nfft)/N; spclkw = 2*abs(spclkw(1:Nfft/2));
spclkp = fft(clkp,Nfft)/N; spclkp = 2*abs(spclkp(1:Nfft/2));
spchpw = fft(chpw,Nfft)/N; spchpw = 2*abs(spchpw(1:Nfft/2));
spchpp = fft(chpp,Nfft)/N; spchpp = 2*abs(spchpp(1:Nfft/2));

figure(1), clf 
subplot(2,1,1), plot(tim,twin,'r')
subplot(2,2,3), plot((ml(1):mh(2))*DF,fwin,'b')
subplot(2,2,4), hold on
plot((ml(1):mh(2))*DF,aw,'k')
plot((ml(1):mh(2))*DF,ap,'m')

figure(2), clf 
subplot(2,2,1), hold on
plot(tim,clkw,'k')
plot(tim,clkp,'r--')
subplot(2,2,2), hold on
plot(tim,chpw,'k')
plot(tim,chpp,'r--')
subplot(2,2,3), hold on
plot(frq,spclkw,'k')
plot(frq,spclkp,'r--')
subplot(2,2,4), hold on
plot(frq,spchpw,'k')
plot(frq,spchpp,'r--')

Dur = 10*1000;
Nfft = 2^nextpow2(Dur/ST);
frq = 1/(ST*2)*linspace(0,1,Nfft/2);

fhp = [0 0.5 1 2 4 8];
filthp = [];
for I = 1:length(fhp)
    filthp = [filthp,(frq>=fhp(I))'];
end
filthp = [filthp;flipud(filthp)];

filtee = 1./sqrt(lcfERB(frq)); filtee = [filtee fliplr(filtee)]';

Flp = 12; flp = lcfNErb2F(lcfF2NErb(Flp)+[-0.5 0.5]);
filtlp = [ones(size(find(frq<=flp(1)))) cos(pi/2*(frq(and(frq>flp(1),frq<=flp(2)))-flp(1))/diff(flp)) zeros(size(find(frq>flp(2))))]; 
filtlp = [filtlp fliplr(filtlp)]';

% White noise now with 12-kHz lowpass filter;
noisew = ifft(fft(randn(Nfft,length(fhp))).*repmat(filtlp,1,length(fhp))); 
noisew = noisew./repmat(sqrt(mean(noisew.^2)),Nfft,1)*10^(-10/20);
noisew = ifft(fft(noisew).*filthp);

% Pink noise now with 12-kHz lowpass filter;
noisep = ifft(fft(randn(Nfft,length(fhp))).*repmat(filtee.*filtlp,1,length(fhp))); 
noisep = noisep./repmat(sqrt(mean(noisep.^2)),Nfft,1)*10^(-10/20);
noisep = ifft(fft(noisep).*filthp);

relLw = zeros(size(fhp)); 
figure(3), clf
c = {'r','g','b','m','c','k'};
subplot(2,1,1), hold on
for I = 1:length(fhp)
    plot(frq,filthp(1:Nfft/2,I).*filtlp(1:Nfft/2),c{I})
end
subplot(2,1,2), hold on
for I = 1:length(fhp)
    jwd = fft(noisew(:,I)); plot(frq,2*abs(jwd(1:Nfft/2)),c{I})
    relLw(I) = 10*log10(sum((2*abs(jwd(1:Nfft/2))).^2)); 
end
% Relative level of hp-filtered white noises relative to unfiltered;
relLw = relLw-relLw(1)

relLp = zeros(size(fhp)); 
figure(4), clf
subplot(2,1,1), hold on
for I = 1:length(fhp)
    plot(frq,filthp(1:Nfft/2,I).*filtee(1:Nfft/2).*filtlp(1:Nfft/2),c{I})
end
subplot(2,1,2), hold on
for I = 1:length(fhp)
    jwd = fft(noisep(:,I)); plot(frq,2*abs(jwd(1:Nfft/2)),c{I})
    relLp(I) = 10*log10(sum((2*abs(jwd(1:Nfft/2))).^2)); 
end
% Relative level of hp-filtered pink noises relative to unfiltered;
relLp = relLp-relLp(1)

% ***** lcfPhase *****
function phase = lcfPhase(f,K,D)

phase = -2*pi*K/(1-D)*f.^(1-D); %Phase delay in radians, taken from Elberling 2007
    
% ***** lcfERB *****
function erb = lcfERB(f)

erb = 24.67*(4.37*f+1)/1000;
    
% ***** lcfF2NErb *****
function nErb = lcfF2NErb(f)

nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);

% ***** lcfNErb2F *****
function f = lcfNErb2F(nErb)

f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
