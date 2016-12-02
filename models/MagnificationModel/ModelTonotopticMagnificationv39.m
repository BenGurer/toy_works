function TonotopicMap = ModelTonotopticMagnification
% Model of Frequency magnification imaged using fMRI
% effect on stimulus spacing on recorded response compared to underlying
% spacing of neuron/voxel characterisic frequency spacing

%% TO DO
%% FUNCTION CREATEFIGURES use if statement and local function to create figures for each stage - input which properties to use

% why does converlution cause the resonse to increase so much max 100 to max 20000
% make sure size of varibles is pre defined - namely each varible in
% TonotopicMap(iTonoMap).sturucture
clc; close all

%% Stimuli Magnification domain functions
StimulusSetSpacing = {str2func('@lcfLinearMagnification'),str2func('@lcfNDLF10'),str2func('@lcfNErb'),str2func('@lcfLogarithmicMagnification')};
StimulusSetFreq = {str2func('@lcfLinMag2Freq'),str2func('@lcfInvNDLF10'),str2func('@lcfInvNErb'),str2func('@lcfLogMag2Freq')};
StimulusSetNames = {'Linear','Frequency Discrimination','Equivalent Rectangular Bandwidth','Logarithmic'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% Stimuli Properties
StimHighFreq = 8;  % highest frequency to test
StimLowFreq = 0.1; % lowest frequency to test
nStimuli = 20;   % number of stimulus
StimdB = 70;    % presentaiton level (dB SPL)

%% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.02; % Lowest frequency of the system
CorticalDistance = 30; % length of gradient in mm
NeuronDensity = 1000; % number of neurons per mm
nNeurons = NeuronDensity .* CorticalDistance; % number of neurons across cortical distance

%% Neuronal population Magnification domain fucntion
% TonotopicMagnification = {str2func('@lcfNDLF'),str2func('@lcfNDLFe'),str2func('@lcfNDLF10'),str2func('@lcfNErb'),str2func('@lcfLogarithmicMagnification')};
% TonoMag2Freq = {str2func('@lcfInvNDLF'),str2func('@lcfInvNDLFe'),str2func('@lcfInvNDLF10'),str2func('@lcfInvNErb'),str2func('@lcfLogMag2Freq')};
% TuningWidth = { str2func('@lcfDLF'),str2func('@lcfDLFe'),str2func('@lcfDLF'),str2func('@lcfErb'),str2func('@lcfErb')};
TonotopicMagnification = {str2func('@lcfNDLF10'),str2func('@lcfNErb'),str2func('@lcfLogarithmicMagnification')};
TonoMag2Freq = {str2func('@lcfInvNDLF10'),str2func('@lcfInvNErb'),str2func('@lcfLogMag2Freq')};
TuningWidth = {str2func('@lcfDLF10'),str2func('@lcfErb'),str2func('@lcfErb')};
% TuningWidth = { str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb')};
TonotopicNames = {'Frequency Discrimination - Log10','Equivalent Rectangular Bandwidth', 'LOG'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

% TonotopicMagnification = {,str2func('@lcfNErb')};
% TonoMag2Freq = {str2func('@lcfInvNErb')};
% % TuningWidth = {str2func('@lcfErb')};
% nr = TWnormalise (1,@lcfErb,@lcfDLF); %(f,TWrFun,TWnFun) % normalise Tuning Width of Tonotopic magnifications
% nr = 1;
% TonotopicMagnificationNames = {'Frequency Discrimination','Equivalent Rectangular Bandwidth'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% System Properites
C = 1; % constant scaling for tuning width - alter Q across all frequencies, for all Tonotopic magnificaiton domains
NoiseLev = 10; % level of noise in system

%% Imaging Properties
% Voxel Population Properties
Resolution = 0.5;  % resolution of imaging in mm
nVoxels = CorticalDistance./Resolution; % number of voxels
VoxelDistance = linspace (1,CorticalDistance,nVoxels);
NeuronsPerVoxel = nNeurons./nVoxels; % number of neurons per voxel
HDspread = 3; % spread of hemodynamic response in mm
%% Model Hemodyanmic Spread
fwhm = 3;% spread of hemodynamic response in mm
sigma = fwhm ./(2 .* (sqrt(2 .* log(2))));
x = -10:(1/NeuronsPerVoxel):10;
Xg = 1 .* exp(-0.5 .* (x./sigma).^2);

%% Figure Properites
% create colours to use for plot each stimulus/voxel property
SpacingColour = {[1 0 0],[0 0.8 0],[0 0 1],[0 0.8 0.8],[0 1 1]};
% use code below to set size of each figure created
% set(f,'OuterPosition',FIG)
FigSize = lcfSetFigSize;
fFontSize = 20;
graphselect = 0; % 1 is for entire frequency range - else zoomed
%% Create indexing varibles
iTonoMap = length(TonotopicMagnification);
iStimSet = length(StimulusSetSpacing);
StimulusSet = zeros (length(StimulusSetSpacing),nStimuli);
StimulusLevel = zeros (length(StimulusSetSpacing),nStimuli);
%% Create stimulus sets
for iStimSet = 1:length(StimulusSetSpacing)
    %% Create Stimulus Set = convert tonotopic magnification units to frequency
    %- linearly space in tonotopic magnification domain between low and high frequency in nStimuli steps
    StimulusSet(iStimSet,:) = StimulusSetFreq{iStimSet}(linspace(StimulusSetSpacing{iStimSet}(StimLowFreq), StimulusSetSpacing{iStimSet}(StimHighFreq), nStimuli));
    StimulusLevel(iStimSet,:) = ones(size(StimulusSet(iStimSet))).*StimdB; % can use to model filter response
end

% TonotopicMap = struct(iTonoMap,1);
%% Present stimulus sets to each tonotopic map
for iTonoMap = 1:length(TonotopicMagnification)
    %% create tonotopic maps
    
    %   Check tuning width normalisation
    nrCheck = func2str(TuningWidth{iTonoMap});
    switch nrCheck
        case 'lcfErb'
            % only normalise DLM p values (filter coeffients) - 3 is the place
            % in the list that ERBs are
            nr = 1;
        otherwise
            nr = TWnormalise (1,@lcfErb,TuningWidth{iTonoMap}); %(f,TWrFun,TWnFun) % normalise Tuning Width of Tonotopic magnifications
    end
    %     if iTonoMap == 4
    %         % only normalise DLM p values (filter coeffients) - 3 is the place
    %         % in the list that ERBs are
    %         nr = 1;
    %     end
    NeuronalProperties = lcfCreateTonotopicProperties(LowFreq,HighFreq,nNeurons,C,TonotopicMagnification{iTonoMap},TonoMag2Freq{iTonoMap},TuningWidth{iTonoMap},nr);    %(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
    TonotopicMap(iTonoMap).CF = NeuronalProperties.CF;
    TonotopicMap(iTonoMap).p = NeuronalProperties.p ;
    TonotopicMap(iTonoMap).Qnr = NeuronalProperties.Qnr ;
    
    for z = 1:nVoxels
        pCF(z) = mean (NeuronalProperties.CF(1 + (NeuronsPerVoxel * (z-1)) : (NeuronsPerVoxel * (z-1)) + NeuronsPerVoxel));
    end
    TonotopicMap(iTonoMap).pCF = pCF;
    
    %% Present each Stimulus Set in turn to neuronal population
    
    NeuronalResponse = zeros(length(StimulusSet),nNeurons);
    NeuronalHemoResponse = zeros(length(StimulusSet),nNeurons);
    VoxelResponse = zeros(length(StimulusSet),nVoxels);
    
    for iStimSet = 1:length(StimulusSetSpacing)
        %         StimFreq = StimulusSet(iStimSet,:);
        %         StimLev = StimulusLevel(iStimSet,:);
        for i = 1:length(StimulusSet)
            NeuronalResponse(i,:) = lcfModelNeuronalResponse(StimulusSet(iStimSet,i),StimulusLevel(iStimSet,i),TonotopicMap(iTonoMap).CF,TonotopicMap(iTonoMap).p ,NoiseLev);
            NeuronalResponse(NeuronalResponse<0) = 0; % remove negative values
            %% Model Hemodyanmic Spread
            NeuronalHemoResponse(i,:) = conv(NeuronalResponse(i,:),Xg,'same');
            for n = 1:nVoxels
                %% fMRI sampling - voxels
                VoxelResponse(i,n) = mean (NeuronalHemoResponse(i,1 + (NeuronsPerVoxel .* (n-1)) : (NeuronsPerVoxel .* (n-1)) + NeuronsPerVoxel));
            end
        end
        %% pCF ESTIMATION
        % weighted average response for each voxel in stimuli domain (stimuli domain because stimuli ID used for calulation)
        ESTpCF = lcfCalculatepCFEstimation (StimulusSet(iStimSet,:),VoxelResponse,StimulusSetSpacing{iStimSet},StimulusSetFreq{iStimSet});
        TonotopicMap(iTonoMap).ESTpCF(iStimSet,:) = ESTpCF';
        
        %% Calculate Root Mean Squared Error
        SE = zeros(1,nVoxels);
        for i = 1:nVoxels
            SE(i)=(pCF(i) - ESTpCF(i)).^2;
        end
        MSE = mean(SE);
        TonotopicMap(iTonoMap).RMSE(iStimSet)= sqrt(MSE);
        
    end
    
%     T = array2table(TonotopicMap(iTonoMap).RMSE,'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'});
%     writetable(T,'RMSE.csv');
end
%% Display results
%% select graphs
if graphselect == 1
    
    % create figure plotting estimated pCF for each magnification domain for
    % each stimulus set with axis limits set by Auditory system min and maximum
    % frequency
    for iTonoMap = 1:length(TonotopicMagnification)
        fTonotopic(iTonoMap) = figure('Color', [1 1 1]);
        set(fTonotopic(iTonoMap),'OuterPosition',FigSize)
        plot (1:nVoxels,TonotopicMap(iTonoMap).pCF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2);
        axis tight, title (TonotopicNames{iTonoMap}, 'FontSize', 20);
        ylim ([LowFreq HighFreq]),xlim ([0 30]);
        hold on
        %     plot ((nNeurons./NeuronsPerVoxel),TonotopicMap(iTonoMap).CF, 'color','b', 'LineWidth', 2);
        for iStimSet = 1:length(StimulusSetNames)
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            %         EstpCF = [zeros(1,xmin), EstpCF(xmin+1:xmax), zeros(1,length(EstpCF) - xmax)];
            scatter (1:nVoxels,EstpCF, 50, SpacingColour{iStimSet}, 'fill');
        end
        set(gca,'YScale','log');
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        fLeg = legend('Actual',(StimulusSetNames{1}),(StimulusSetNames{2}),(StimulusSetNames{3}),(StimulusSetNames{4}),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
    
    
else
    % create figure plotting estimated pCF for each magnification domain for
    % each stimulus set with axis limits set by stimulus set min and maximum
    for iTonoMap = 1:length(TonotopicMagnification)
        fTonotopicZoomPrune(iTonoMap) = figure('Color', [1 1 1]);
        set(fTonotopicZoomPrune(iTonoMap),'OuterPosition',FigSize)
        plot (VoxelDistance,TonotopicMap(iTonoMap).pCF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2);
        axis tight, title (TonotopicNames{iTonoMap}, 'FontSize', 20);
        xmin = find(TonotopicMap(iTonoMap).pCF >= StimLowFreq,1);
        xmax = find(TonotopicMap(iTonoMap).pCF >= StimHighFreq,1);  % highest frequency to test
                 xminDis = (xmin/ nVoxels) * max(VoxelDistance);
         xmaxDis = (xmax/ nVoxels) * max(VoxelDistance);
        TonotopicMap(iTonoMap).VoxelsWithinStimLim(iTonoMap) = xmax - xmin; % need to convert from voxels to mm+
        TonotopicMap(iTonoMap).CorticalDistanceWithinStimLim(iTonoMap) = xmaxDis - xminDis; 
        ylim ([0.2 8]), xlim ([xminDis xmaxDis]);
        hold on
        for iStimSet = 1:length(StimulusSetNames)
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            EstpCF = [zeros(1,xmin), EstpCF(xmin+1:xmax), zeros(1,length(EstpCF) - xmax)];
            scatter (VoxelDistance,EstpCF, 50, SpacingColour{iStimSet}, 'fill');
        end
        set(gca,'YScale','log');
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        fLeg = legend('Actual',(StimulusSetNames{1}),(StimulusSetNames{2}),(StimulusSetNames{3}),(StimulusSetNames{4}),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
end

%% Calculate Root Mean Squared Error for voxels within stimulus frequency limits

% data.RMSEStimLim = zeros(1,(iTonoMap));
for iTonoMap = 1:length(TonotopicMagnification)
    MSE = zeros(1,(iStimSet));
    % prune by pCF highest and lowest frequency in stimulus sets
    pCF = TonotopicMap(iTonoMap).pCF;
    xmin = find(pCF >= StimLowFreq,1);
    xmax = find(pCF >= StimHighFreq,1);  % highest frequency to test
    pCF = pCF(xmin+1:xmax);
    TonotopicMap(iTonoMap).VoxelsWithinStimLim = xmax -xmin;
    for iStimSet = 1:length(StimulusSetNames)
        EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
        EstpCF = EstpCF(xmin+1:xmax);
        SE = zeros(iStimSet,length(pCF));
        for i = 1:length(pCF)
            SE(iStimSet,i)= (pCF(i) - EstpCF(i)).^2;
        end
        MSE(iStimSet) = mean(SE(iStimSet,:));
        TonotopicMap(iTonoMap).RMSEPruned(iStimSet) = sqrt(MSE(iStimSet));
    end
end
% T = array2table(TonotopicMap(iTonoMap).RMSEPruned(iStimSet),'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'});
% writetable(T,'RMSEStimPruned.csv');
%% Stimulus sets plot
fStim = figure('Color', [1 1 1]);
% create figure to plot stimulus properties
for iStimSet = 1:length(StimulusSetSpacing)
    for i = 1:length((StimulusSet(iStimSet,:)))
        figure (fStim)% plot to stimulus properties figure
        subplot(2,2,iStimSet); % create and index subplot
        plot(repmat(StimulusSet(iStimSet,i),1,2),[min(ylim) StimulusLevel(iStimSet,i)], 'color',[(StimulusSet(iStimSet,i))./max(StimulusSet(iStimSet,:)) 0.5 0.5], 'LineWidth', 10); % at each stimulus frequency create lines the level of each stimulus presenation level
        % set plot properties
        set(gca,'XScale','log');
        title (StimulusSetNames{iStimSet}, 'FontSize', fFontSize), xlabel('Frequency (kHz)', 'FontSize', fFontSize), ylabel('Presentation Level (dB SPL)', 'FontSize', fFontSize);
        ylim ([0 StimulusLevel(iStimSet,i)]), xlim ([LowFreq 20]);
        hold on
    end
end

%% plot Charaterisic Frequencies for each magnification scale
fCF = figure('Color', [1 1 1]);
for iTonoMap = 1:length(TonotopicMagnification)
    plot (1:nNeurons,TonotopicMap(iTonoMap).CF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2)
    hold on
    set(gca,'YScale','log');
    xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('CF (kHz)', 'FontSize', fFontSize);
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    fLeg = legend((TonotopicNames),'Location','SouthEast');
    set(fLeg, 'FontSize', fFontSize./2);
end

%% Plot tuning widths as a function of CF
fTuningvsCF = figure('Color', [1 1 1]);
for iTonoMap = 1:length(TonotopicMagnification)
    plot (TonotopicMap(iTonoMap).Qnr ,TonotopicMap(iTonoMap).CF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2)
    hold on
    set(gca,'YScale','log','XScale','log');
    xlabel('Tuning Width (Q)', 'FontSize', fFontSize), ylabel('CF (kHz)', 'FontSize', fFontSize);
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    fLeg = legend((TonotopicNames),'Location','SouthEast');
    set(fLeg, 'FontSize', fFontSize./2);
end
end
%% Local Functions
function properties = lcfCreateTonotopicProperties(LowFreq,HighFreq,nNeurons,C,TonotopicMagnification,TonoMag2Freq,TuningWidth,nr)
%(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
TonotopicFrequencySpacing = linspace(TonotopicMagnification(LowFreq),TonotopicMagnification(HighFreq),nNeurons); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
% Characterisic Frequencies
properties.CF = TonoMag2Freq(TonotopicFrequencySpacing);   % neuron Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
% Tuning Width
Q = TuningWidth(properties.CF); % Neuron frequency tuning width for each Voxel Charactersic Frequency (for each neuron)(Q)
properties.Qnr = Q ./ nr;    % normalise Tuning width at 1 kHz by ERB
properties.p = C*4*properties.CF./properties.Qnr; % Filter coefficient of each Voxel % VTuningWidth - population tuning width / populationERB
% properties.por = C*4*properties.CF./Q; % filter coefficient not normalised
end
function Response = lcfModelNeuronalResponse (StimulusSet,StimulusLevel,CF,p,NoiseLev)
StimInt = 10.^(StimulusLevel/10);
Int = zeros(1,length(CF)); % Response intensity for each Voxel
% loop through all neurons
for I = 1:length(CF)
    for II = 1:length(StimulusSet)
        g = abs((StimulusSet(II)-CF(I))/CF(I));
        fw=(1+(p(I)*g)).*exp(-p(I)*g); %Two parameter roex
        Int(I) = Int(I)+fw*StimInt(II);
    end
end
% NoiseLev = 10;
noise = randn (1,length(CF)).*NoiseLev;
Response = 10*log10(max(Int,1)); % Threshold intensity response and convert to dB
Response = Response + noise; %% set negative values to 0
% Response = Int; % no threshold
% Response = max(Int,1);
end
function ESTpCF = lcfCalculatepCFEstimation (StimFreq, VoxRes,StimMagDomain,StimMagDomain2Freq)
ESTpCF = zeros(length(VoxRes),1);
StimFreqStimDomain = StimMagDomain(StimFreq);
for i = 1:length(VoxRes)
    ESTpCF (i) = sum((StimFreqStimDomain'.*VoxRes(:,i)))./sum(VoxRes(:,i)); % need to convert to linear spacing - convert to mag domain then convert back for characteristic frequencies
end
ESTpCF(isnan(ESTpCF)) = 0 ;

ESTpCF = StimMagDomain2Freq(ESTpCF);
%% need to convert back to frequency domain
end
function x = lcfLinearMagnification(f)
x = f;
end
function f = lcfLinMag2Freq(x)
f = x;
end
function X = lcfLogarithmicMagnification(f)
X = log10(f);
end
function f = lcfLogMag2Freq(x)
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


function DLF = lcfDLF10 (f)
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = 10^(a*sqrt(f)+B))
% for 40 dB SL
a = 0.026;
b = -0.533;
f = f.*1000;
DLF = 10.^(a*sqrt(f)+b);
DLF = DLF./1000;
end
function nDLF = lcfNDLF10 (f)
f = f.*1000;    % convert from kHz to Hz
% log df = exp(a*sqrt(f)+B)) Weir 1977
% for 40 dB SL
a = 0.026;
b = -0.533;
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = 10^(a*sqrt(f)+B)) = df/dx
% dx/df = 1/(10^(a*sqrt(x)+b))
% integral = x
% x = -(2^(1-b-a sqrt(x)) 5^(-b-a sqrt(x)) (1+a sqrt(x) log(10)))/(a^2 log^2(10))
nDLF = -(2.^(1-b-a .* sqrt(f)) .* 5.^(-b-a .* sqrt(f)) .* (1+a .* sqrt(f) .* log(10)))./(a.^2 .*(log(2)+log(5)).^2);
nDLF = nDLF./1000;
end
function f = lcfInvNDLF10 (nDLF)
fs = 0:0.01:20;
fs = fs.*1000;  % convert from kHz to Hz
a = 0.026;
b = -0.533;
nDLF = nDLF.*1000;
DLF = -(2.^(1-b-a .* sqrt(fs)) .* 5.^(-b-a .* sqrt(fs)) .* (1+a .* sqrt(fs) .* log(10)))./(a.^2 .*(log(2)+log(5)).^2);
f = interp1(DLF,fs,nDLF,'spline');
f = f./1000; % convert from Hz to kHz
end


function DLF = lcfDLFe (f)
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = e^(a*sqrt(f)+B))
% for 40 dB SL
a = 0.026;
b = -0.533;
f = f.*1000;
DLF = exp(a*sqrt(f)+b);
DLF = DLF./1000;

end
function nDLF = lcfNDLFe (f)
f = f.*1000;    % convert from kHz to Hz
% log df = exp(a*sqrt(f)+B)) Weir 1977
% for 40 dB SL
a = 0.026;
b = -0.533;
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = e^(a*sqrt(f)+B)) = df/dx
% dx/df = % 1/e^(a*sqrt(x)+b)
% integral = x
% x = -(2 e^(-b-a sqrt(x)) (1+a sqrt(x)))/a^2
nDLF = -((2.*(a .* sqrt(f) + 1) .* exp(a .* (-sqrt(f)) - b))./ a.^2);
nDLF = nDLF./1000;
end
function f = lcfInvNDLFe (nDLF)
fs = 0:0.01:20;
fs = fs.*1000;  % convert from kHz to Hz
a = 0.026;
b = -0.533;
nDLF = nDLF.*1000;
DLF = -((2.*(a .* sqrt(fs) + 1) .* exp(a .* (-sqrt(fs)) - b))./ a.^2);
f = interp1(DLF,fs,nDLF,'spline');
f = f./1000; % convert from Hz to kHz
end


function nr = TWnormalise (f,TWrFun,TWnFun)
% nr = normalising ratio
% f = frequenct to reference to (kHz)
% TWrFun = function which defines the reference tuning width
% TWnFun = function which defines the tuning width to normalise

r = TWrFun(f);
n = TWnFun(f);

nr = n ./ r;

end
function FigSize = lcfSetFigSize
% Get monitor size to make figures full screen
f = figure ('Visible','off',...
    'Menu', 'None', ...
    'Color', [0.5, 0.5, 0.5]);

set(0,'Units','pixels') ;
p = get(0, 'MonitorPositions');
scnsize = p(1,:);
position = get(f,'Position');
outerpos = get(f,'OuterPosition');
borders = outerpos - position;
edge = -borders(1)/2;

close (f)

FigSize = [scnsize(1),...
    scnsize(2),...
    scnsize(3) - edge,...
    scnsize(4)];
end
