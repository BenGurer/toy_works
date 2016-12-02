function data = LogVsLinModel
figure(1), clf

StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10;
% compF = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps)
compF = linspace(StimLowFreq,StimHighFreq,StimSteps)
% compF = 1;
compLev= ones(size(compF));
C=0.06;
plotC = {'k','b','r','g','y'}; % Cell array of colros.

plotC = repmat(plotC,1,round(length(compF))./(StimSteps./2));

for i = 1:length(compF)
    excPat  = funROEX(compF(i),compLev(i),C,@funF2NErb,@funNErb2F,@funErb);
    figure(1)
    plot(excPat.VoxelFreqC,excPat.excLev,'Color',plotC{i}), set(gca,'XScale','log'), axis tight, xlabel('fc (kHz)'), ylabel('excLev (dB)')
    hold on
    line(repmat(compF(i),1,2),[min(ylim) compLev(i)],'Color',plotC{i})
end
end

function VoxelResponse = funROEX(compF,compLev,C,magDomain,magScale2F,magFilterWidth)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

HighFreq = 20;
LowFreq = 0.05;

Gradient = 30; % length of gradient in mm
Resolution = 1;  % resolution of image in mm
nVoxels = Gradient/Resolution;
VoxelResponse.magScale = linspace(magDomain(LowFreq),magDomain(HighFreq),nVoxels); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
% excPat.nErb = linspace(funF2NErb(LowFreq),funF2NErb(HighFreq),nVoxels);
% excPat.nErb = funF2NErb(LowFreq):NeuronSpacing:funF2NErb(HighFreq);
VoxelResponse.VoxelFreqC = magScale2F(VoxelResponse.magScale);   % Voxel Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
% excPat.fc = funNErb2F(excPat.magScale);
VoxelResponse.erb = magFilterWidth(VoxelResponse.VoxelFreqC);
% excPat.erb = funErb(excPat.VoxelFreqC);
VoxelResponse.p = C*4*VoxelResponse.VoxelFreqC./VoxelResponse.erb;

VoxelResponse.excInt = zeros(1,length(VoxelResponse.VoxelFreqC));
for I = 1:length(VoxelResponse.VoxelFreqC)
    for II = 1:length(compF)
        g = abs((compF(II)-VoxelResponse.VoxelFreqC(I))/VoxelResponse.VoxelFreqC(I));
        fw=(1+(VoxelResponse.p(I)*g)).*exp(-VoxelResponse.p(I)*g); %Two parameter roex
        VoxelResponse.excInt(I) = VoxelResponse.excInt(I)+fw*compInt(II);
    end
end
VoxelResponse.excLev = 10*log10(max(VoxelResponse.excInt,1)); % Threshold intensity response and convert to dB
end
function erb = funErb(f)
% Erb (in kHz) according to Glasberg and Moore 1990 as a function of frequency, f (in kHz);
% f can be a vector.

erb = 24.67*(4.37*f+1)/1000;
end
function f = funNErb2F(nErb)
% Converts Erb number (nErb) to frequency, f (in kHz); nErb can be a vector.

f = (10.^(nErb*24.67*4.37/(1000*log(10)))-1)/4.37;
end
function nErb = funF2NErb(f)
% Converts frequency, f (in kHz), to Erb number (nErb) according to Glasberg and Moore 1990;
% f can be a vector; nErb = integral from f = 0 to f of 1/erb(f).

nErb = 1000*log(10)/(24.67*4.37)*log10(4.37*f+1);
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