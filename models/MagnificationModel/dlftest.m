function dlftest
% clc; close all;
lf = 0.02;
hf = 20;
% r = 20;
% f = linspace(lf,hf,r);
step = 0.01;
f = lf:step:hf;
n =length(f)
nDLF = lcfNDLF10 (f);
freturn = lcfInvNDLF10 (nDLF);

% fminus = freturn - f;
% figure
% plot(fminus)
% figure
% plot (f,freturn)

stepcomp = (hf - lf) /n;
dlf = lcfDLF10 (f);
dlfdif = diff(lcfNDLF10 (f))./step;
dlfrec = 1./dlfdif;
% ndlfoverf = (lcfNDLF10 (f))./f;

figure
% plot (f,dlf)
% hold on
plot (f(1:end-1)+step./2,dlfrec, 'r')
    set(gca,'YScale','log','XScale','log');

end

function DLF = lcfDLF10 (f)
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = 10^(a*sqrt(f)+B))
% for 40 dB SL
a = 0.026;
b = -0.533;
f = f.*1000;
DLF = 10.^(a*sqrt(f)+b);
DLF = DLF./1000;
end
function nDLF = lcfNDLF10 (f)
f = f.*1000;    % convert from kHz to Hz
% log df = exp(a*sqrt(f)+B)) Weir 1977
% for 40 dB SL
a = 0.026;
b = -0.533;
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = 10^(a*sqrt(f)+B)) = df/dx
% dx/df = 1/(10^(a*sqrt(x)+b))
% integral = x
% x = -(2^(1-b-a sqrt(x)) 5^(-b-a sqrt(x)) (1+a sqrt(x) log(10)))/(a^2 log^2(10))
nDLF = -(2.^(1-b-a .* sqrt(f)) .* 5.^(-b-a .* sqrt(f)) .* (1+a .* sqrt(f) .* log(10)))./(a.^2 .*(log(2)+log(5)).^2);
% nDLF = nDLF./1000;
end
function f = lcfInvNDLF10 (nDLF)
fs = 0:0.001:25;
fs = fs.*1000;  % convert from kHz to Hz
a = 0.026;
b = -0.533;

DLF = -(2.^(1-b-a .* sqrt(fs)) .* 5.^(-b-a .* sqrt(fs)) .* (1+a .* sqrt(fs) .* log(10)))./(a.^2 .*(log(2)+log(5)).^2);
f = interp1(DLF,fs,nDLF,'spline');
f = f./1000; % convert from Hz to kHz
end