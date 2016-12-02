function VMD = fMRIFreqMagModel
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
% create fucntions for each mag domian
% ceate loop to index magdomians there fucntion handle
% plot neuron response as line, voxel as point

%% Neuronal population Magnification domain fucntion
% VoxelMagSpace = {str2func('@lcfNErb'),str2func('@FreqSelectMagDom')};
% VoxelMagSpace2Freq = {str2func('@lcfInvNErb'),str2func('@FreqSelectMagDom2Fre')};
% VoxelFilterWidth = {str2func('@lcfErb'),str2func('@FreqSelectMagDom2Freq')};

%% Neuronal population Magnification domain fucntion
VoxelMagSpace = {str2func('@lcfNErb')};
VoxelMagSpace2Freq = {str2func('@lcfInvNErb')};
VoxelFilterWidth = {str2func('@lcfErb')};

%% Stimuli Magnification domain functions
% magDomain,magScale2F,magFilterWidth
StimMagDomain = {str2func('@LinMagDom'),str2func('@LogMagDom'),str2func('@lcfNErb'),str2func('@FreqSelectMagDom')};
StimMagDomain2Freq = {str2func('@LinMagDom2Freq'),str2func('@LogMagDom2Freq'),str2func('@lcfInvNErb'),str2func('@FreqSelectMagDom2Freq')};

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

VMD = cell(length(VoxelMagSpace),1);
SMD = struct('VoxRes', ones(StimNum,nVoxels), ...
    'VCharaFreq', ones(nVoxels), ...
    'fMRIdata', ones(nVoxels));
Res = cell(length(StimMagDomain),1);
for iVMD = 1:length(VoxelMagSpace)
    % Spacing Properites
    VoxelProp = VoxelProperties(LowFreq,HighFreq,nVoxels,C,VoxelMagSpace{iVMD},VoxelMagSpace2Freq{iVMD},VoxelFilterWidth{iVMD});    %     (LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
    
    for iSMD = 1:length(StimMagDomain)
        %% create stimuli frequencies for each spacing
        StimFreq = StimMagDomain2Freq{iSMD}(linspace(StimMagDomain{iSMD}(StimLowFreq),StimMagDomain{iSMD}(StimHighFreq),StimNum)); % frequencies of stimuli - 3D matrix for band limited noise
        
        %% create stimuli gain structure
        StimLev = ones(size(StimFreq)); % can use to model filter response
        
        %% loop for each voxel
        
        for i = 1:length(StimFreq)
            SMD.VoxRes(i,:) = VoxelResponse (StimFreq(i),StimLev,VoxelProp.VFreqC,VoxelProp.p);
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
    VMD{iVMD} = Res;
end
%% plot Voxel characterisic frequencies - actual vs measured
figure
subplot(2,2,1);
plot (VMD{1,1}{1,1}.VCharaFreq)
subplot(2,2,2);
plot (VMD{1,1}{3,1}.VCharaFreq)
subplot(2,2,3);
plot (VMD{1,1}{4,1}.VCharaFreq)
subplot(2,2,4);
plot (VMD{1,1}{3,1}.VCharaFreq)
% figure
% subplot(2,2,1);
% scatter (1:nVoxels,VoxelProp.VFreqC, 'r'), axis tight
% hold on
% scatter (1:nVoxels,data.MeasuredVoxFreqC{1,1}, 'b','MarkerFaceColor','b'),set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Characterisic Frequency (kHz)'),title ('Log'),legend('Actual','Measured','Location','SouthEast')% ,legend('Actual','Measured','Location','SouthEastOutside')
% % strStimFreq= ['Stimuli = ',num2str(ymin)];
% % text(xmin,ymin,strStimFreq,'HorizontalAlignment','left');
% 
% subplot(2,2,2);
% scatter (1:nVoxels,VoxelProp.VFreqC, 'r'), axis tight
% hold on
% scatter (1:nVoxels,data.MeasuredVoxFreqC{1,2}, 'b','MarkerFaceColor','b'),set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Characterisic Frequency (kHz)'),title ('Lin'),legend('Actual','Measured','Location','SouthEast')%,legend('Actual','Measured','Location','SouthEastOutside')
% 
% subplot(2,2,3);
% scatter (1:nVoxels,VoxelProp.VFreqC, 'r'), axis tight
% hold on
% scatter (1:nVoxels,data.MeasuredVoxFreqC{1,3}, 'b','MarkerFaceColor','b'),set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Characterisic Frequency (kHz)'),title ('ERB'),legend('Actual','Measured','Location','SouthEast')
% 
% data.VoxelCharFreq = VoxelProp.VFreqC;
% data.StimFreq = StimFreq;

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
function lin = LinMagDom(f)
lin = f;
end
function f = LinMagDom2Freq(lin)
f = lin;
end
function log = LogMagDom(f)
log = f;
end
function f = LogMagDom2Freq(log)
f = log;
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
function FreqSelTunWid = FreqSelTuningWidth (f)
FreqSelTunWid = f;
end
function FreqSel = FreqSelectMagDom(f)
FreqSel = f;
end
function f = FreqSelectMagDom2Freq(FreqSel)
f = FreqSel;
end