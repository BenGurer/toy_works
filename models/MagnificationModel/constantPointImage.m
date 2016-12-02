function constantPointImage
close all
f = 0.02:0.1:20;
Da = 0;
Db = 30;
Fa = 2;
Fb = 8;

for i = 1:2
    if i == 1
  G = @lcfNErb;
g = @lcfErb;
Gi = @lcfInvNErb;
    elseif i == 2
          G = @lcfNDLF10;
g = @lcfDLF10;
Gi = @lcfInvNDLF10;
    end
    
B = ((Db - Da) * (G(Fb))) / ((G(Fa)) - (G(Fb)));
A = (Da - B)/G(Fa);
D = A * (G(f)) + B;

% figure
% plot(f,D)

TWn = 0.1;
Fn = 2;

C = (TWn * Fn) / g(Fn);
Bc = -(C * (G(Fn)));
Dc = C .* (G(f)) + Bc;
Fc = Gi(((Dc) - Bc)./C);

erbc = C * g(f);
figure
plot(f,Dc)
figure
plot(Dc,Fc)
figure
plot(f,erbc)
end
end
function erb = lcfErb(f)
% ***** lcfErb *****
% ERBs as per Glasberg and Moore (1990);
A = 24.7/1000; B = 4.37;
erb = A*(B*f+1);
end
function nerb = lcfNErb(f)
% ***** lcfNErb *****
% Converts frequency to ERB number;
A = 24.7/1000; B = 4.37;
nerb = 1/(A*B)*log(B*f+1);
end
function f = lcfInvNErb(nerb)
% ***** lcfInvNErb *****
% Converts ERB number to frequency;
A = 24.7/1000; B = 4.37;
f = 1/B*(exp(A*B*nerb)-1);
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
end
function f = lcfInvNDLF10 (nDLF)
%% interpolates to estimate frequency as a function of number of DLFs
fs = 0:0.01:20;
fs = fs.*1000;  % convert from kHz to Hz
a = 0.026;
b = -0.533;
DLF = -(2.^(1-b-a .* sqrt(fs)) .* 5.^(-b-a .* sqrt(fs)) .* (1+a .* sqrt(fs) .* log(10)))./(a.^2 .*(log(2)+log(5)).^2);
f = interp1(DLF,fs,nDLF,'spline');
f = f./1000; % convert from Hz to kHz
end
