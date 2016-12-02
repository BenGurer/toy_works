function f = funInvNDLF (nDLF)
close all
%% define frequency resolution
fs = 0.02:0.1:20;
% convert from kHz to Hz
fs = fs.*1000;
cortexDistance = 30;
% AuditoryLowFreq = AuditoryLowFreq .*1000;
% constants for frequency discrimination
% SL is dB SL
a = 0.0214; k = -0.15; m = 5.056; SL = 35;
% y = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% equation from Nelson
% DLM = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% integral for 1/exp(a .* sqrt(f) + k + m.*(SL.^-1))
% describes how distance changes as a function of characteristic frequency
% Fc
% df = -((2.*(a.*sqrt(f) + 1) .* exp (-(a.* SL .* sqrt(f) + k .*SL + m)./SL))./a.^2);

% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = exp(a*sqrt(f)+B))
 % for 40 dB SL
a = 0.026;
b = -0.533;
fd = exp(a*sqrt(fs)+b);
df = 1./(exp(a*sqrt(fs)+b));
% integral of 1/e^(a*sqrt(f)+b)
% d = ((-((2.*(a.*sqrt(fs) +1).*exp(a.*(-sqrt(fs))-b)./a.^2))+3000)./2645).*30;
% 1/e^(a*sqrt(x)+b)
% integral
% -(2 e^(-b-a sqrt(x)) (1+a sqrt(x)))/a^2
d = -((2.*(a .* sqrt(fs) + 1) .* exp(a .* (-sqrt(fs)) - b))./ a.^2);
% 1/(10^(a*sqrt(x)+b))
% integral
% -(2^(1-b-a sqrt(x)) 5^(-b-a sqrt(x)) (1+a sqrt(x) log(10)))/(a^2 log^2(10))
dlog10 = -(2.^(1-b-a .* sqrt(fs)) .* 5.^(-b-a .* sqrt(fs)) .* (1+a .* sqrt(fs) .* log(10)))./(a.^2 .* (log(10)).^2);
d = (((d - min(d)))./(max(d)-min(d))).*cortexDistance;
dlog10 = (((dlog10 - min(dlog10)))./(max(dlog10)-min(dlog10))).*cortexDistance;
% inverse f = 1/(a.^2 .* log(10).6^2)

A = 24.7/1000; B = 4.37;
% A = 24.7; B = 4.37;
nerb = 1/(A*B)*log(B*fs+1);
% nerb = (log(A.*B.*fs +A))./A.*B;
% nerb = (log(A.*(B.*fs + 1)))./A.*B;
nerb = (((nerb - min(nerb)))./(max(nerb)-min(nerb))).*cortexDistance;

lind = fs;
lind = (((lind - min(lind)))./(max(lind)-min(lind))).*cortexDistance;

logd = log(fs);
logd = (((logd - min(logd)))./(max(logd)-min(logd))).*cortexDistance;

% d = -((2.*(a.*sqrt(fs) +1).*exp(a.*(-sqrt(fs)-b))./a.^2));
% d = (((d - min(d)))./(max(d)-min(d))).*cortexDistance;

dq = 0:0.01:30; % distance of cortex
% x = cortical distance = d
% v = characterisic frequency = f
fq = interp1(d,fs,dq,'spline');
% fq = interp1(f,d,dq,'spline');
% figure
% semilogx (f,df)
% axis tight
f = interp1(d,fs,nDLF,'spline')
f = f./1000
figure
subplot (2,2,1)
plot (fs,fd)
subplot (2,2,2)
plot (fs,df)
subplot (2,2,3)
plot (fs,d)
hold on
plot (fs,logd,'r')
plot (fs,lind,'b')
plot (fs,nerb,'g')
plot (fs,dlog10, 'p')

subplot (2,2,4)
% plot (dq,fq,'o')
plot(d,fs, 'o',dq,fq)
hold on
% figure
plot (logd,fs,'r')
% hold on
plot (lind,fs,'b')
plot (nerb,fs,'g')
plot (dlog10,fs, 'p')
% semilogy (dq,fq)
% semilogx (d,f, 'o',dq,fq)
% semilogx (fq,dq)
% semilogx (f,d, 'o',dq,fq)

% (1+ProductLog(1/2 a^2 e^(-1+k+m/S) y))^2/a^2
% dq = 0:1:30;
% fq = ((1+lambertw(1/2 .* a.^2 .* exp(-1+k+m./SL) .* d)).^2) ./ a.^2;
end
