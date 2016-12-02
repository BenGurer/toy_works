function r = TWnormalise (f)
clc; close all
% f = 1;

fDLF = f .*1000 % convert to Hz to work with formula
a = 0.026;
b = -0.533;
DLF = 10.^(a*sqrt(f)+b);

% ***** lcfErb *****
% ERBs as per Glasberg and Moore (1990);
A = 24.7/1000; B = 4.37;
erb = A*(B*f+1);
r = DLF ./ erb;

normDLF = DLF ./ r;
end

function nr = TWnormalise (f,TWrFun,TWnFun)
% nr = normalising ratio
% f = frequenct to reference to (kHz)
% TWrFun = function which defines the reference tuning width
% TWnFun = function which defines the tuning width to normalise

r = TWrFun(f);
n = TWnFun(f);

nr = n ./ r;

end
function DLF = lcfDLF (f)
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = 10^(a*sqrt(f)+B))
 % for 40 dB SL
a = 0.026;
b = -0.533;
DLF = 10.^(a*sqrt(f)+b);

end
function erb = lcfErb(f)
% ***** lcfErb *****
% ERBs as per Glasberg and Moore (1990);
A = 24.7/1000; B = 4.37;
erb = A*(B*f+1);
end