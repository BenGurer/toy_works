function nErb = funF2NErbGM90(f)
% Converts frequency, f (in kHz), to Erb number (nErb) according to Glasberg and Moore 1990; 
% f can be a vector; nErb = integral from f = 0 to f of 1/erb(f).

nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);
