function data = ModelTonotopticMagnification
% Model of Frequency magnification imaged using fMRI
% effect on stimulus spacing on recorded response compared to underlying
% spacing of neuron/voxel characterisic frequency spacing

%% TO DO
clc; close all

%% Stimuli Magnification domain functions
StimulusSetSpacing = {str2func('@LinMagDom'),str2func('@lcfNDLF'),str2func('@lcfNErb'),str2func('@LogMagDom')};
StimulusSetFreq = {str2func('@LinMagDom2Freq'),str2func('@lcfInvNDLF'),str2func('@lcfInvNErb'),str2func('@LogMagDom2Freq')};
StimulusSetNames = {'Linear','Frequency Discrimination','Equivalent Rectangular Bandwidth','Logarithmic'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% Stimuli Properties
StimHighFreq = 8;  % highest frequency to test
StimLowFreq = 0.25; % lowest frequency to test
nStimuli = 10;   % number of stimulus
StimdB = 70;    % presentaiton level (dB SPL)


%% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.02; % Lowest frequency of the system
CorticalDistance = 30; % length of gradient in mm
NeuronDensity = 10000; % number of neurons per mm
nNeurons = NeuronDensity .* CorticalDistance; % number of neurons across cortical distance

%% Neuronal population Magnification domain fucntion
TonotopicMagnification = {str2func('@lcfNDLF'),str2func('@lcfNErb')};
TonoMag2Freq = {str2func('@lcfInvNDLF'),str2func('@lcfInvNErb')};
TuningWidth = { str2func('@lcfDLF'),str2func('@lcfErb')};
nr = TWnormalise (1,@lcfErb,@lcfDLF); %(f,TWrFun,TWnFun) % normalise Tuning Width of Tonotopic magnifications
TonotopicMagnificationNames = {'Frequency Discrimination','Equivalent Rectangular Bandwidth'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% System Properites
C = 1; % constant scaling for tuning width - alter Q across all frequencies, for all Tonotopic magnificaiton domains
NoiseLev = 10; % level of noise in system

%% Imaging Properties
% Voxel Population Properties
Resolution = 1;  % resolution of imaging in mm
nVoxels = CorticalDistance./Resolution; % number of voxels
NeuronsPerVoxel = nNeurons./nVoxels; % number of neurons per voxel
HDspread = 3; % spread of hemodynamic response in mm

%% Figure Properites
% create colours to use for plot each stimulus/voxel property
SpacingColour = {[1 0 0],[0 0.8 0],[0 0 1],[0 0.8 0.8]};
% use code below to set size of each figure created
% set(f,'OuterPosition',FIG)
FigSize = lcfSetFigSize;
fFontSize = 20;

%% Create indexing varibles
iTonoMap = length(TonotopicMagnification);
iStimSet = length(StimulusSetSpacing);

%% Create stimulus sets
for iStimSet = 1:length(StimulusSetSpacing)
    
    %% Create Stimulus Set = convert tonotopic magnification units to frequency
    %- linearly space in tonotopic magnification domain between low and high frequency in nStimuli steps
    StimulusSet(iStimSet) = StimulusSetFreq{iStimSet} (linspace (StimulusSetSpacing{iStimSet}(StimLowFreq), StimulusSetSpacing{iStimSet}(StimHighFreq), nStimuli));
    StimulusLevel(iStimSet) = ones(size(StimulusSet(iStimSet))).*StimdB; % can use to model filter response
    
end

%% Present stimulus sets to each tonotopic map

for iTonoMap = 1:length(TonotopicMagnification)
    %% create tonotopic maps
    NeuronalProperties = CreateTonotopicProperties(LowFreq,HighFreq,nNeurons,C,TonotopicMagnification{iTonoMap},TonoMag2Freq{iTonoMap},TuningWidth{iTonoMap},nr);    %(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)      
    TonotopicMap(iTonoMap).CF = NeuronalProperties(iTonoMap);
    for iStimSet = 1:length(StimulusSetSpacing)
        for i = 1:length(StimulusSet)
            NeuronalResponse(i,:) = modelNeuronalResponse(StimulusSet(i),StimulusLevel,NeuronalProperties.CF,NeuronalProperties.p,NoiseLev);
            NeuronalResponse(NeuronalResponse<0) = 0; % remove negative values
            
            %% FUNCTION CREATEFIGURES use if statement and local function to create figures for each stage - input which properties to use
        end
        %% fMRI sampling - voxels
        for n = 1:nVoxels
            VoxelResponse(n) = mean (NeuronalResponse(1 + (NeuronsPerVoxel .* (n-1)) : (NeuronsPerVoxel .* (n-1)) + NeuronsPerVoxel));
        end
        
        %% Model Hemodyanmic Spread
        % gaussian filter with 3 mm spread - 1.5 either side
%         VoxelResponse

        %% pCF ESTIMATION
        % weighted average response for each voxel in stimuli domain (stimuli domain because stimuli ID used for calulation)
        pCF = pCFEstimation (StimulusSet{iStimSet},VoxelResponse,StimulusSetSpacing{iStimSet},StimulusSetFreq{iStimSet});
     TonotopicMap(iTonoMap).(iStimSet).pCF
    end
end


%% Local Functions
function neuron = CreateTonotopicProperties(LowFreq,HighFreq,nNeurons,C,TonotopicMagnification,TonoMag2Freq,TuningWidth,nr);    %(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
    
neuron.TonotopicMagnification = linspace(TonotopicMagnification(LowFreq),TonotopicMagnification(HighFreq),nNeurons); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
% Characterisic Frequencies
neuron.CF = TonoMag2Freq(TonotopicMagnification);   % neuron Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
% Tuning Width
neuron.Q = TuningWidth(CF); % Neuron frequency tuning width for each Voxel Charactersic Frequency (for each neuron)(Q)
neuron.Q = neuron.Q ./ nr;
neuron.p = C*4*neuron.CF./neuron.Q; % Filter coefficient of each Voxel % VTuningWidth - population tuning width / populationERB
end

end