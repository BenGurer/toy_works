function ERB = getERB(f)
%
%   usage: getERB(f)
%      by: Ben Gurer
%    date: 23/01/2017
% purpose: calculate ERB for given values of frequency (kHz)
% 
% discription:
% ERBs as per Glasberg and Moore (1990). Frequency input as kHz.
ERB = zeros(1,length(f));

A = 24.7; B = 4.37;
ERB = A*(B*f+1);
end