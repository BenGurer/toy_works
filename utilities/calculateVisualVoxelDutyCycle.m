%% Hansen 2004
% First figure out eccentricity of Hansen's stimuli

r1 = 3;
r2 = 19;

% S = cos(C*r^(1/6))
C = (9*pi)/(r2^(1/6)-r1^(1/6))

a1 = C*r1^(1/6)
a2 = C*r2^(1/6)

a = linspace(a1,a2,10)

r = (a./C).^6

% Apply non-linear scaling to the stimuli used by Dumoulin to determine the
% centre of the stimulus range in non-linear space (where it is most likely
% his estiamte of tuning width is correct)
duR1 = 0.25
duR2 = 14


duA1 = C*duR1^(1/6)
duA2 = C*duR2^(1/6)
aCentre =  (duA2 + duA1)/2

dur = (aCentre./C).^6

% the smallest stimulus used by Hansen is more eccentric than the the
% location that Dumoulins estimate is most likely correct

% one stimulus = 3 degrees width
sigma = 1;
fwhm = 2.355* sigma % convert sigma to FWHM
