function data = fMRIFreqMagModel
% Model of Frequency magnification imaged using fMRI
% effect on stimulus spacing on recorded response compared to underlying
% spacing of neuron/voxel characterisic frequency spacing

%% Stimulus Properties
StimHighFreq = 12;
StimLowFreq = 0.25;
StimSteps = 10;

StimFreq = cell(1,3);
StimFreq {1,1} = logspace(log10(StimLowFreq),log10(StimHighFreq),StimSteps);
StimFreq {1,2}= linspace(StimLowFreq,StimHighFreq,StimSteps);    % frequencies of stimuli - 3D matrix for band limited noise
StimFreq {1,3}= linspace(lcfNErb(StimLowFreq),lcfNErb(StimHighFreq),StimSteps);
StimFreq {1,3}= lcfInvNErb(StimFreq{1,3});
StimLev = cell(1,3);
for n = 1:3
StimLev {1,n} = ones(size(StimFreq{1,n})); % Level of stimulus frequencies - can use to filter
end

%% Voxel Population Properties
% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.05; % Lowest frequency of the system
Gradient = 30; % length of gradient in mm
% Imaging Properties
Resolution = 1;  % resolution of imaging in mm
nVoxels = Gradient/Resolution; % number of voxels

% Spacing Properites
C = 0.06; % constant scaling for voxel tuning width
VoxelProp = VoxelProperties(LowFreq,HighFreq,nVoxels,C,@lcfNErb,@lcfInvNErb,@lcfErb);

%% fMRI response data
fMRIdata = cell(1,3);
data.MeasuredVoxFreqC = cell(1,3);
for n = 1:3
fMRIdata {1,n} = zeros(length(StimFreq{1,n}),length(VoxelProp.VmagScale)); % create matrix to fill with loop - x =number of stimulus y = number of voxels
for i = 1:length(StimFreq{1,n})
fMRIdata{1,n}(i,:) = HDresponse(StimFreq{1,n}(i),StimLev{1,n}(i),VoxelProp.VFreqC,VoxelProp.p);
%fMRI data = StimFreq(10),VoxPos(30)
end
%% weight average response for each voxel in stimuli domain (stimuli domain because stimuli ID used for calulation)
data.MeasuredVoxFreqC {1,n}= zeros(nVoxels,1);
for i = 1:nVoxels
data.MeasuredVoxFreqC{1,n}(i) = sum((StimFreq{1,n}'.*fMRIdata{1,n}(:,i)))./sum(fMRIdata{1,n}(:,i));
end
data.MeasuredVoxFreqC{1,n}(isnan(data.MeasuredVoxFreqC{1,n})) = 0 ;
end
%% plot Voxel characterisic frequencies
figure
subplot(2,1,1); 
plot (1:nVoxels,VoxelProp.VFreqC, 'r'), axis tight
hold on
plot (data.MeasuredVoxFreqC, 'b','LineWidth',2),set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Characterisic Frequency (kHz)')
hleg1 = legend('Actual','Measured','Location','SouthEastOutside');

subplot(2,1,2); 
scatter (1:nVoxels,VoxelProp.VFreqC, 'r'), axis tight
hold on
scatter (1:nVoxels,data.MeasuredVoxFreqC, 'b','MarkerFaceColor','b'),set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Characterisic Frequency (kHz)'),legend('Actual','Measured','Location','SouthEastOutside')

% plot (fMRIdata)

end
function Response = HDresponse (StimF,StimLev,VoxelFreqC,p)
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

end
function VoxelProp = VoxelProperties(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
VoxelProp.VmagScale = linspace(magDomain(LowFreq),magDomain(HighFreq),nVoxels); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
VoxelProp.VFreqC = magScale2F(VoxelProp.VmagScale);   % Voxel Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
VoxelProp.VTuningWidth = magFilterWidth(VoxelProp.VFreqC); % Voxel frequency tuning width for each Voxel Charactersic Frequency (for each voxel)(Bandwidth in?)
VoxelProp.p = C*4*VoxelProp.VFreqC./VoxelProp.VTuningWidth; % Filter coefficient of each Voxel
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