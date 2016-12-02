function data = fMRIFreqMagModel
% Model of Frequency magnification imaged using fMRI
% effect on stimulus spacing on recorded response compared to underlying
% spacing of neuron/voxel characterisic frequency spacing

%% TO DO
% Add spread of HMDR - apply to voxel response
% Add noise
% Add frequency selectivity spacing
% create bandpass noise stimuli
% each voxel contain many tuning curves - loop current code for each voxel
% and take weighted average - then convovle with gaussian Hymodynmaic
% spread response
% plot superimposed stimuli responses
% plot neuron response as line, voxel as point
% 30 random numbers - create coiefient 

clc; close all

%% Neuronal population Magnification domain fucntion
% VoxelMagSpace = {str2func('@lcfNErb'),str2func('@FreqDisMagDom')};
% VoxelMagSpace2Freq = {str2func('@lcfInvNErb'),str2func('@FreqDisMagDom2Fre')};
% VoxelFilterWidth = {str2func('@lcfErb'),str2func('@FreqDisMagDom2Freq')};

%% Neuronal population Magnification domain fucntion
VoxelMagSpace = {str2func('@lcfNErb')};
VoxelMagSpace2Freq = {str2func('@lcfInvNErb')};
VoxelFilterWidth = {str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb')};

%% Stimuli Magnification domain functions
% magDomain,magScale2F,magFilterWidth
StimMagDomain = {str2func('@LinMagDom'),str2func('@LogMagDom'),str2func('@lcfNErb'),str2func('@FreqDisMagDom')};
StimMagDomain2Freq = {str2func('@LinMagDom2Freq'),str2func('@LogMagDom2Freq'),str2func('@lcfInvNErb'),str2func('@FreqDisMagDom2Freq')};
StimSpacingNames = {'Lin','Log','ERB','FreqDis'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% Stimuli Properties
StimHighFreq = 12;  % highest frequency to test
StimLowFreq = 0.25; % lowest frequency to test
StimNum = 10;   % number of stimulus
%% Voxel Population Properties
% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.05; % Lowest frequency of the system
Gradient = 30; % length of gradient in mm
nNeurons = 100000; % number of neurons per mm
% Imaging Properties
Resolution = 1;  % resolution of imaging in mm
nVoxels = Gradient/Resolution; % number of voxels
NeuronDensity = nNeurons .* Resolution; % number of neurons in each voxel
HDspread = 3; % spread of hemodynamic response in mm

% Spacing Properites
C = 0.06; % constant scaling for voxel tuning width

VoxelMagSpace = StimMagDomain ;
VoxelMagSpace2Freq = StimMagDomain2Freq;

VMD = cell(length(VoxelMagSpace),1);
SMD = struct('VoxRes', ones(StimNum,nVoxels), ...
    'VCharaFreq', ones(nVoxels), ...
    'fMRIdata', ones(nVoxels));
Res = cell(length(StimMagDomain),1);
PlotIndex = [0 4 8 12]

for iVMD = 1:length(VoxelMagSpace)
    % Spacing Properites
    data.VoxelProp(iVMD)  = VoxelProperties(LowFreq,HighFreq,nVoxels,C,VoxelMagSpace{iVMD},VoxelMagSpace2Freq{iVMD},VoxelFilterWidth{iVMD});    %     (LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
    
    for iSMD = 1:length(StimMagDomain)
        %% create stimuli frequencies for each spacing
        StimFreq = StimMagDomain2Freq{iSMD}(linspace(StimMagDomain{iSMD}(StimLowFreq),StimMagDomain{iSMD}(StimHighFreq),StimNum)); % frequencies of stimuli - 3D matrix for band limited noise
        
        %% create stimuli gain structure
        StimLev = ones(size(StimFreq)); % can use to model filter response
        
        %% loop for each voxel
         figure (1)
        for i = 1:length(StimFreq)
            SMD.VoxRes(i,:) = VoxelResponse (StimFreq(i),StimLev,data.VoxelProp(iVMD).VFreqC,data.VoxelProp(iVMD).p);
            subplot(length(VoxelMagSpace),length(StimMagDomain),iSMD + PlotIndex(iVMD));
            plot(data.VoxelProp(iVMD).VFreqC,SMD.VoxRes(i,:)), set(gca,'XScale','log'),title (StimSpacingNames{iSMD}), xlabel('Actual Voxel Characterisic Frequency (kHz)'), ylabel('Response Level'), axis tight
            hold on
            line(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)])
        end
        
        % weight average response for each voxel in stimuli domain (stimuli domain because stimuli ID used for calulation)
        %         for i = 1:nVoxels
        SMD.VCharaFreq = VoxelCharacterFrequency (StimFreq, SMD.VoxRes,StimMagDomain{iSMD},StimMagDomain2Freq{iSMD});
        %         end
        
        %% fMRI response data
        %% convolve with Hemodynamic spread
        %         SMD.fMRIdata(iSMD) = fMRIResponse (SMD.VCharaFreq(iSMD),HDspread);
        
        Res{iSMD}= SMD;
        
    end
    data.VMD{iVMD,1} = Res;
end  
    %% plot Voxel characterisic frequencies - actual vs measured
    for iVMD = 1:length(VoxelMagSpace)
        figure (2)
        for iSMD = 1:length(StimMagDomain)
            subplot(length(VoxelMagSpace),length(StimMagDomain),iSMD + PlotIndex(iVMD));
            plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'r'), axis tight
            hold on
            scatter (1:nVoxels,data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq)
            set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Estimated Characterisic Frequency (kHz)'),title (StimSpacingNames{iSMD}),legend('Actual','Measured','Location','SouthEast')
        end
    end
    
% end

end
function Response = VoxelResponse (StimF,StimLev,VoxelFreqC,p)
StimInt = 10.^(StimLev/10);
Int = zeros(1,length(VoxelFreqC)); % Response intensity for each Voxel
for I = 1:length(VoxelFreqC)
    for II = 1:length(StimF)
        g = abs((StimF(II)-VoxelFreqC(I))/VoxelFreqC(I));
        fw=(1+(p(I)*g)).*exp(-p(I)*g); %Two parameter roex
        Int(I) = Int(I)+fw*StimInt(II);
    end
end
NoiseLev = 0.1;
noise = randn (1,length(VoxelFreqC)).*NoiseLev;
Int = Int + noise;
Response = 10*log10(max(Int,1)); % Threshold intensity response and convert to dB
% Response = Int; % no threshold
% Response = max(Int,1);
end
function VCharaFreq = VoxelCharacterFrequency (StimFreq, VoxRes,StimMagDomain,StimMagDomain2Freq)
VCharaFreq = zeros(length(VoxRes),1);
StimFreqStimDomain = StimMagDomain(StimFreq);
for i = 1:length(VoxRes)
    VCharaFreq (i) = sum((StimFreqStimDomain'.*VoxRes(:,i)))./sum(VoxRes(:,i)); % need to convert to linear spacing - convert to mag domain then convert back for characteristic frequencies
end
VCharaFreq(isnan(VCharaFreq)) = 0 ;

VCharaFreq = StimMagDomain2Freq(VCharaFreq);
%% need to convert back to frequency domain
end
function Response = fMRIResponse (StimF,StimLev,VoxelFreqC,p)
% hemodynamic spread
%% Parkes (2005) at 3 T
%% Find paper for 7T - around 2.3 mm
% FWHM GE = 3.9 +/- 0.7 mm, SE = 3.4 +/- 0.8
StimInt = 10.^(StimLev/10);
Int = zeros(1,length(VoxelFreqC)); % Response intensity for each Voxel
for I = 1:length(VoxelFreqC)
    for II = 1:length(StimF)
        g = abs((StimF(II)-VoxelFreqC(I))/VoxelFreqC(I));
        fw=(1+(p(I)*g)).*exp(-p(I)*g); %Two parameter roex
        Int(I) = Int(I)+fw*StimInt(II);
    end
end
Response = 10*log10(max(Int,1)); % Threshold intensity response and convert to dB
% Response = Int; % no threshold
% Response = max(Int,1);
end
function VoxelProp = VoxelProperties(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
VoxelProp.VmagScale = linspace(magDomain(LowFreq),magDomain(HighFreq),nVoxels); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
VoxelProp.VFreqC = magScale2F(VoxelProp.VmagScale);   % Voxel Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
VoxelProp.VTuningWidth = magFilterWidth(VoxelProp.VFreqC); % Voxel frequency tuning width for each Voxel Charactersic Frequency (for each voxel)(Bandwidth in?)
VoxelProp.p = C*4*VoxelProp.VFreqC./VoxelProp.VTuningWidth; % Filter coefficient of each Voxel
end
function y = LinMagDom(f)
y = f;
end
function f = LinMagDom2Freq(x)
f = x;
end
function y = LogMagDom(f)
y = log10(f);
end
function f = LogMagDom2Freq(x)
f = 10.^(x);
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
function DLF= FreqDisTuningWidth (f)

%% log DLF = a'sqrt(F) + k' + m'(SL^-1) % Nelson et al (1983)
% Nelson et al k' = -0.15, m' = 5.056, a' = 0.0214

%% Weir (1977) = log DLF = a sqrt(F) + b % depended on dB SL

a = 0.0214; k = -0.15; m = 5.056; SL = 40;

DLF = log10(a(sqrt(f)) + k + m(SL^-1));


% FreqSelTunWid = f;
end
function y = FreqDisMagDom(f)
%% number of Freq Selectivity units
a = 0.0214; k = -0.15; m = 5.056; SL = 35;
y = exp(a .* sqrt(f) + k + m.*(SL.^-1));
end
function f = FreqDisMagDom2Freq(x)
a = 0.0214; k = -0.15; m = 5.056; SL = 35;
f = ((log(x) - (k + m.*(SL.^-1)))./a).^2;
end