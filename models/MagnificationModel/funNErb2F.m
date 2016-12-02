function f = funNErb2FGM90(nErb)
% Converts Erb number (nErb) to frequency, f (in kHz); nErb can be a vector.

f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
