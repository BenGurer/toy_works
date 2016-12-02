function Response = LogVsLinModel
clf
f = figure(1);
        axes1 = axes('Parent',f,'YGrid','on',...
            'XScale','linear',...
            'XMinorTick','on',...
            'XMinorGrid','on',...
            'XGrid','on');
StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10;
% compF = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps)
compF = linspace(StimLowFreq,StimHighFreq,StimSteps);
% compF = 1;
compLev= ones(size(compF));
C=0.06;
plotC = {'k','b','r','g','y'}; % Cell array of colros.

plotC = repmat(plotC,1,round(length(compF))./(StimSteps./2));
Response = cell(size(compF));
for i = 1:length(compF)
    Response{i}  = funROEX(compF(i),compLev(i),C,@funF2NErb,@funNErb2F,@funErb);
    
%     plot(axes1, excPat.VoxelFreqC,excPat.excLev,'Color',plotC{i}), set(gca,'XScale','log'), axis tight, xlabel('Frequency Gradient (mm)'), ylabel('Nom Response (arb)')
%         plot(axes1, Response{i}.VoxelFreqC,Response{i}.excLev,'Color',plotC{i}), axis tight, xlabel('Frequency Gradient (mm)'), ylabel('Nom Response (arb)')
%     hold on
%     line(repmat(compF(i),1,2),[min(ylim) compLev(i)],'Color',plotC{i})
end

plot (compF,Response {:,1}.excLev(6))

end

function Response = funROEX(compF,compLev,C,magDomain,magScale2F,magFilterWidth)
% compF/compLev are the frequencyies and levels of the stimulus frequency
% components; Ctip/Ctail are teh widths of the tip and tail roex filters; W
% is the relative gain of the tail filter;

compInt = 10.^(compLev/10);

HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.05; % Lowest frequency of the system

Gradient = 30; % length of gradient in mm
Resolution = 1;  % resolution of image in mm
nVoxels = Gradient/Resolution;
Response.magScale = linspace(magDomain(LowFreq),magDomain(HighFreq),nVoxels); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
Response.VoxelFreqC = magScale2F(Response.magScale);   % Voxel Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
Response.VoxelTuningWidth = magFilterWidth(Response.VoxelFreqC);
Response.p = C*4*Response.VoxelFreqC./Response.VoxelTuningWidth;

Response.excInt = zeros(1,length(Response.VoxelFreqC));
for I = 1:length(Response.VoxelFreqC)
    for II = 1:length(compF)
        g = abs((compF(II)-Response.VoxelFreqC(I))/Response.VoxelFreqC(I));
        fw=(1+(Response.p(I)*g)).*exp(-Response.p(I)*g); %Two parameter roex
        Response.excInt(I) = Response.excInt(I)+fw*compInt(II);
    end
end
Response.excLev = 10*log10(max(Response.excInt,1)); % Threshold intensity response and convert to dB
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