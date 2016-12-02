function erb = funErbGM90(f)
% Erb (in kHz) according to Glasberg and Moore 1990 as a function of frequency, f (in kHz); 
% f can be a vector.

erb = 24.67*(4.37*f+1)/1000;
