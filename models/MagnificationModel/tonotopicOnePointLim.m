function tonotopicOnePointLim
clc; close all;
lf = 0.02;
hf = 20;
step = 0.01;

f = lf:step:hf;

df = lcfErb(f);

figure
semilogx(f,df,'c')
hold on
index = find(f==2);
nr = (10* f(index)/100) / df(index);
figure
dfn = df * nr;
semilogx(f,dfn,'g')
figure
dxdf = 1 ./ dfn;
semilogx(f,dxdf,'r')

x = cumsum(dxdf);

figure
semilogx(f,x,'b')

xn = x - x(index);

figure
semilogx(f,xn)

figure
nerb = lcfNErb(f);
semilogx(f,nerb,'c')
hold on
nerbn = nerb ./nr;
semilogx(f,nerbn)
nerbnoffset = nerbn - nerbn(index);

figure
semilogx(f,nerbnoffset)
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