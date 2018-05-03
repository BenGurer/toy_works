function [percent NoisePercent] = calTuningWidthandNoise
Fn = 2;
dataNERB = linspace(lcfNErb(0.250),lcfNErb(8),200);
dataFrequencies = lcfInvNErb(dataNERB);
diffERB = mean(diff(dataNERB));
fun = @(x,mu,sigma) 1 * exp(-(x - mu).^2/2/sigma^2);
mu = lcfNErb(Fn);
sig = [13.5];    % sigma in erb space
sig = mean(sig)*diffERB;
q = integral(@(x)fun(x,mu,sig),-100,100); % erb in erb space
erb = lcfInvNErb(mu + q/2) - lcfInvNErb(mu - q/2);
% percent = (erb .* 100)/Fn;
percent = (erb/Fn) .* 100;

MaxRes = [6];
NoiseRes = [1.5322];
NoisePercent = (NoiseRes ./ MaxRes).*100;
NoisePercent = mean(NoisePercent);
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