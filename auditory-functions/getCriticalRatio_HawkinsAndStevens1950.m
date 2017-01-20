function [f_measured_hz, CriticalRatio] = getCriticalRatio_HawkinsAndStevens1950
%
%   usage: getCriticalRatio_HawkinsAndStevens1950
%      by: Ben Gurer
%    date: 19/01/2017
% purpose: output the critical ratio values as determined by Hawkins And Stevens (1950)
% 
% discription:
% values taken from Hawkins And Stevens (1950) The masking of pure tones and of speech by white noise (fig. 6)
% "Ratio between the monaural masked threshold of a pure tone and the level per cycle of the masking noise measured at the frequency of the pure tone".

 
f_measured_hz = [100, 125, 175, 250, 350, 500, 700, 1000, 1400, 2000, 2800, 4000, 5600, 7000, 8000, 9000];
CriticalRatio = [19, 17.75, 17.5, 16.25, 16.5 17.25, 17.75, 18.5, 19.25, 20.5, 22.5, 25, 25.5, 26.75, 27, 28.5];

% 
% [0.125 0.25 0.5 1 1.5 2 3 4 6 8];
% [17.75 16.3 17.25 18.5 19.25 20.5 22.5 25.1 26 27];
% figure; semilogx(f_measured_hz,CriticalRatio)