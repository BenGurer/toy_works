function nDLF = funNDLF (f)
%% define frequency resolution
% fs = 0:0.1:20;
% convert from kHz to Hz
f = f.*1000;
% constants for frequency discrimination
% SL is dB SL
% a = 0.0214; k = -0.15; m = 5.056; SL = 35;
% y = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% equation from Nelson
% DLM = exp(a .* sqrt(f) + k + m.*(SL.^-1));
% integral for 1/exp(a .* sqrt(f) + k + m.*(SL.^-1))
% describes how distance changes as a function of characteristic frequency
% Fc
% df = -((2.*(a.*sqrt(f) + 1) .* exp (-(a.* SL .* sqrt(f) + k .*SL + m)./SL))./a.^2);

% log df = exp(a*sqrt(f)+B)) Weir 1977

 % for 40 dB SL
a = 0.026;
b = -0.533;
% fd = exp(a*sqrt(fs)+b);
% df = 1./(exp(a*sqrt(fs)+b));
% integral of 1/e^(a*sqrt(f)+b)
% d = ((-((2.*(a.*sqrt(fs) +1).*exp(a.*(-sqrt(fs)-b))./a.^2))+3000)./2645).*30;

d = -((2.*(a.*sqrt(f) +1).*exp(a.*(-sqrt(f)-b))./a.^2));
nDLF = -((2.*(a.*sqrt(f) +1).*exp(a.*(-sqrt(f)-b))./a.^2));
% dq = 0:00.1:30; % distance of cortex
% x = cortical distance = d
% v = characterisic frequency = f
% fq = interp1(d,fs,dq,'spline');
% fq = interp1(f,d,dq,'spline');
% figure
% semilogx (f,df)
% axis tight
% nDLF = interp1(d,fs,f,'spline')
% figure
% subplot (2,2,1)
% plot (fs,fd)
% subplot (2,2,2)
% plot (fs,df)
% subplot (2,2,3)
% plot (fs,d)
% subplot (2,2,4)
% plot (dq,fq)
% semilogy (dq,fq)
% semilogx (d,f, 'o',dq,fq)
% semilogx (fq,dq)
% semilogx (f,d, 'o',dq,fq)

% (1+ProductLog(1/2 a^2 e^(-1+k+m/S) y))^2/a^2
% dq = 0:1:30;
% fq = ((1+lambertw(1/2 .* a.^2 .* exp(-1+k+m./SL) .* d)).^2) ./ a.^2;
end
