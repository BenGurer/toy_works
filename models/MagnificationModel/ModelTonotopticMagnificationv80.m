function TonotopicMap = ModelTonotopticMagnification
% Model of Frequency magnification imaged using fMRI
% effect on stimulus spacing on recorded response compared to underlying
% spacing of neuron/voxel characterisic frequency spacing

%% TO DO
%% FUNCTION CREATEFIGURES use if statement and local function to create figures for each stage - input which properties to use

% why does converlution cause the resonse to increase so much max 100 to max 20000
% make sure size of varibles is pre defined - namely each varible in - TonotopicMap(iTonoMap).sturucture
% Bandlimit RMSE
% calculate RMSE in tonotopic magnifiction domain
% plot graphs for each section to check

% ability to diferenciaite two frequencies isn't the same as tuning width
% of neruons - response from many neurons - more neurons that response the
% more accuritate - ERB based on measurements in notched noise so only it
% is reflective of tuning width

% point image per mm needed to voxel based fMRI sampling - allow to be set different for
% each magnification - massively different scales

% number of point images per mm determined by percentage

% point image needs to be whole numbers evenly spaced to work with model

% sample x distance from 2 kHz
% vary number of Point images per mm

% round to whole values and set step size for point image

% sample percentage of point images

% convert SPL to SL for Frequency Discrimination equation

% list of conditions to model

% stimulus set

% auditory properties

clc; close all

%% Stimuli Magnification domain functions
StimulusSetSpacing = {str2func('@lcfNDLF10'),str2func('@lcfNErb'),str2func('@lcfLogarithmicMagnification'),str2func('@lcfLinearMagnification')};
StimulusSetFreq = {str2func('@lcfInvNDLF10'),str2func('@lcfInvNErb'),str2func('@lcfLogMag2Freq'),str2func('@lcfLinMag2Freq')};
StimulusSetNames = {'Frequency Discrimination','Equivalent Rectangular Bandwidth','Logarithmic','Linear'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% Stimuli Properties
StimHighFreq = 8;  % highest frequency to test
StimLowFreq = 0.25; % lowest frequency to test
nStimuli = 10;   % number of stimulus
StimdB = 70;    % presentaiton level (dB SPL)

%% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.02; % Lowest frequency of the system
FreqStep = 0.001; % resolution of auditory system
f = LowFreq:FreqStep:HighFreq; % Frequencies representated by auditory system
Fn = 2;     % Frequency to calibrate tonotopic maps (kHz)
TWn = 0.1; % Tuning width at Fn in perceptage
CorticalDistance = 30; % length of gradient in mm
NeuronDensity = 1000; % number of neurons per mm
nNeurons = NeuronDensity .* CorticalDistance; % number of neurons across cortical distance
NeuronDistance = linspace (0,CorticalDistance,nNeurons);

%% Neuronal population Magnification domain fucntion
TonotopicMagnification = {str2func('@lcfNDLF10'),str2func('@lcfNErb')};
TonoMag2Freq = {str2func('@lcfInvNDLF10'),str2func('@lcfInvNErb')};
TuningWidth = {str2func('@lcfDLF10'),str2func('@lcfErb')};


% TuningWidth = {str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb')};
TonotopicNames = {'Frequency Discrimination','Equivalent Rectangular Bandwidth', 'Logarithmic'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables

%% System Properites
% C = 0.01; % constant scaling for tuning width - alter Q across all frequencies, for all Tonotopic magnificaiton domains
C = 1;
NoiseLev = 10; % level of noise in system

%% Imaging Properties
% Voxel Population Properties
Resolution = 0.5;  % resolution of imaging in mm
nVoxels = CorticalDistance/Resolution; % number of voxels
VoxelDistance = linspace (0,CorticalDistance,nVoxels);
NeuronsPerVoxel = nNeurons./nVoxels; % number of neurons per voxel
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
rangelow = 0.5;
rangehigh = 4;

%% Select graphs
% graphselect = 0; % 1 is for entire frequency range - else zoomed
fStimulusSel = 1;
fCFSel = 1;
fTuningWidthSel = 1;
fNeuronResSel = 1;
fNeuronHemoResSel = 1;
fVoxelResSel = 1;
fESTpCFSel = 1;
fESTpCFPrunedSel = 0;
fRMSErange = 0;
fStimulusvsTonotopicSel = 1;

%% Create indexing varibles - not used
iTonoMap = length(TonotopicMagnification);
iStimSet = length(StimulusSetSpacing);

%% Create stimulus sets
%% Predefine varibles to save stimulus sets
StimulusSet = zeros (length(StimulusSetSpacing),nStimuli);
StimulusLevel = zeros (length(StimulusSetSpacing),nStimuli);
for iStimSet = 1:length(StimulusSetSpacing)
    %% Create Stimulus Set = convert tonotopic magnification units to frequency
    %- linearly space in tonotopic magnification domain between low and high frequency in nStimuli steps
    StimulusSet(iStimSet,:) = StimulusSetFreq{iStimSet}(linspace(StimulusSetSpacing{iStimSet}(StimLowFreq), StimulusSetSpacing{iStimSet}(StimHighFreq), nStimuli));
    StimulusLevel(iStimSet,:) = ones(size(StimulusSet(iStimSet))).*StimdB; % can use to model filter response
end

%% Present stimulus sets to each tonotopic map

%% Predefine varibles to store model data
TonotopicMap = struct('SpacingOrder',{'DLM','ERB','Log'},...
    'CF',zeros(1,nNeurons),...
    'p',zeros(1,nNeurons),...
    'df',zeros(1,nNeurons),...
    'pCF',zeros(length(StimulusSetSpacing),nVoxels),...
    'ESTpCF',zeros(length(StimulusSetSpacing),nVoxels),...
    'RMSE',zeros(1,length(StimulusSetSpacing)),...
    'RMSEPruned',zeros(1,length(StimulusSetSpacing)),...
    'VoxelsWithinStimLim',[],...
    'CorticalDistanceWithinStimLim',[]);

for iTonoMap = 1:length(TonotopicMagnification)
    %% create tonotopic maps
    %% Define neuronal properties and save to output variable
    NeuronalProperties = lcfCreateTonotopicProperties(Fn,TWn,f,nNeurons,TuningWidth{iTonoMap},TonotopicMagnification{iTonoMap},TonoMag2Freq{iTonoMap},C);
    TonotopicMap(iTonoMap).D = NeuronalProperties.D;
    TonotopicMap(iTonoMap).CF = NeuronalProperties.CF;
    TonotopicMap(iTonoMap).p = NeuronalProperties.p ;
    TonotopicMap(iTonoMap).df = NeuronalProperties.df ;
    
    %% fMRI
    % sample constant number of point images

    for z = 1:nVoxels
        pCF(z) = mean (NeuronalProperties.CF(1 + (NeuronsPerVoxel * (z-1)) : (NeuronsPerVoxel * (z-1)) + NeuronsPerVoxel));
    end
    TonotopicMap(iTonoMap).pCF = pCF;
    
    
    %% Present each Stimulus Set in turn to neuronal population
    % Predefine variables
    NeuronalResponse = zeros(length(StimulusSet),length(TonotopicMap(iTonoMap).CF));
        NeuronalHemoResponse = zeros(length(StimulusSet),length(TonotopicMap(iTonoMap).CF));
        VoxelResponse = zeros(length(StimulusSet),nVoxels);
    for iStimSet = 1:length(StimulusSetSpacing)
        for i = 1:length(StimulusSet)
            %% Present each stimulus from stimulus set in turn to neuronal population
            % Output response of all neurons to one stimulus presenation
            NeuronalResponse(i,:) = lcfModelNeuronalResponse(StimulusSet(iStimSet,i),StimulusLevel(iStimSet,i),TonotopicMap(iTonoMap).CF,TonotopicMap(iTonoMap).p ,NoiseLev);
            NeuronalResponse(NeuronalResponse<0) = 0; % remove negative values
                                    NoisePercent = 10/100;
                        NoiseLev = StimdB*NoisePercent;
noise = randn (1,length(NeuronalResponse(i,:))).*NoiseLev;
NeuronalResponse(i,:) = NeuronalResponse(i,:) + noise; 
NeuronalResponse(i,:) = max(0,NeuronalResponse(i,:));
            %% Model Hemodyanmic Spread
            % convolve response with guassian function - sigma is set to FWHM of hemodymanic spread
                        NeuronalHemoResponse(i,:) = conv(NeuronalResponse(i,:),Xg,'same');

            for n = 1:nVoxels
                %% fMRI sampling - voxels
                % Assuming neurons are equaly spaced - sample neuronal response by imaging resolution to give voxel response
                VoxelResponse(i,n) = mean (NeuronalHemoResponse(i,1 + (NeuronsPerVoxel .* (n-1)) : (NeuronsPerVoxel .* (n-1)) + NeuronsPerVoxel));
            end
               NoisePercent = 10/100;
                        NoiseLev = (10*10^4)*NoisePercent;
noise = randn (1,length(VoxelResponse(i,:))).*NoiseLev;
VoxelResponse(i,:) = VoxelResponse(i,:) + noise; 
VoxelResponse(i,:) = max(0,VoxelResponse(i,:));             
        end
        
        TonotopicMap(iTonoMap).Response(iStimSet).NeuronalResponse(:,:) = NeuronalResponse;
                TonotopicMap(iTonoMap).Response(iStimSet).NeuronalHemoResponse(:,:) = NeuronalHemoResponse;
        TonotopicMap(iTonoMap).Response(iStimSet).VoxelResponse(:,:) = VoxelResponse;
        %         %% CF ESTIMATION
        %         % weighted average response for each voxel in stimuli domain (stimuli domain because stimuli ID used for calulation)
        ESTpCF = lcfCalculatepCFEstimation (StimulusSet(iStimSet,:),VoxelResponse,StimulusSetSpacing{iStimSet},StimulusSetFreq{iStimSet});
        TonotopicMap(iTonoMap).ESTpCF(iStimSet,:) = ESTpCF';
        
        %         %% Calculate Root Mean Squared Error for all voxels
        SE=(log10(ESTpCF') - log10(pCF)).^2;
        MSE = mean(SE);
        TonotopicMap(iTonoMap).RMSE(iStimSet)= sqrt(MSE);
    end
    
end
%% Display results
%% select graphs
%% Stimulus sets plot
if fStimulusSel == 1
    fStim = figure('Color', [1 1 1]);
    set(fStim,'OuterPosition',FigSize)
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
end
%% plot Charaterisic Frequencies for each magnification scale
if fCFSel == 1
    fCF = figure('Color', [1 1 1]);
    set(fCF,'OuterPosition',FigSize)
    for iTonoMap = 1:length(TonotopicMagnification)
        subplot (length(TonotopicMagnification), 1, iTonoMap)
        plot (TonotopicMap(iTonoMap).D,TonotopicMap(iTonoMap).CF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2)
        %         hold on
        set(gca,'YScale','log');
        xlabel('Point Image', 'FontSize', fFontSize), ylabel('CF (kHz)', 'FontSize', fFontSize);
        ylim ([LowFreq HighFreq]), xlim ([min(TonotopicMap(iTonoMap).D) max(TonotopicMap(iTonoMap).D)]);
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        fLeg = legend((TonotopicNames{iTonoMap}),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
    fCF = figure('Color', [1 1 1]);
    set(fCF,'OuterPosition',FigSize)
    for iTonoMap = 1:length(TonotopicMagnification)
        subplot (length(TonotopicMagnification), 1, iTonoMap)
        plot (TonotopicMap(iTonoMap).CF,TonotopicMap(iTonoMap).D, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2)
        %         hold on
        set(gca,'XScale','log');
        ylabel('Point Image', 'FontSize', fFontSize), xlabel('CF (kHz)', 'FontSize', fFontSize);
        xlim ([LowFreq HighFreq])
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        fLeg = legend((TonotopicNames{iTonoMap}),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
end

if fTuningWidthSel == 1
    %% Plot tuning widths as a function of CF
    fTuningvsCF = figure('Color', [1 1 1]);
    set(fTuningvsCF,'OuterPosition',FigSize)
    for iTonoMap = 1:length(TonotopicMagnification)
        plot (TonotopicMap(iTonoMap).CF,TonotopicMap(iTonoMap).df, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2)
        hold on
        set(gca,'YScale','log','XScale','log');
        ylabel('Tuning Width (kHz)', 'FontSize', fFontSize), xlabel('CF (kHz)', 'FontSize', fFontSize);
        xlim ([LowFreq HighFreq])
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        fLeg = legend((TonotopicNames),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
        
    end
end

if fNeuronResSel == 1;
    %% plot Neuronal response
    for iTonoMap = 1:length(TonotopicMagnification)
        fNeuronRes(iTonoMap) = figure('Color', [1 1 1]);
        set(fNeuronRes(iTonoMap),'OuterPosition',FigSize)
        for iStimSet = 1:length(StimulusSetNames)
            subplot(2,2,iStimSet);
            for i = 1 : length(StimulusSet)
                area(TonotopicMap(iTonoMap).D,TonotopicMap(iTonoMap).Response(iStimSet).NeuronalResponse(i,:), 'FaceColor',[(StimulusSet(iStimSet,i)./max(StimulusSet(iStimSet,:))) 0.5 0.5])
                hold on
                
            end
        end
    end
end

if fNeuronHemoResSel == 1;
    %% plot voxel response
    for iTonoMap = 1:length(TonotopicMagnification)
        fNeuronHemoRes(iTonoMap) = figure('Color', [1 1 1]);
        set(fNeuronHemoRes(iTonoMap),'OuterPosition',FigSize)
        for iStimSet = 1:length(StimulusSetNames)
            subplot(2,2,iStimSet);
            for i = 1 : length(StimulusSet)
                area(TonotopicMap(iTonoMap).D,TonotopicMap(iTonoMap).Response(iStimSet).NeuronalHemoResponse(i,:), 'FaceColor',[(StimulusSet(iStimSet,i)./max(StimulusSet(iStimSet,:))) 0.5 0.5])
                hold on

            end
        end
    end
end

if fVoxelResSel == 1;
    %% plot voxel response
    for iTonoMap = 1:length(TonotopicMagnification)
        fVoxelRes(iTonoMap) = figure('Color', [1 1 1]);
        set(fVoxelRes(iTonoMap),'OuterPosition',FigSize)
        for iStimSet = 1:length(StimulusSetNames)
            subplot(2,2,iStimSet);
            for i = 1 : length(StimulusSet)
                area(1:nVoxels,TonotopicMap(iTonoMap).Response(iStimSet).VoxelResponse(i,:), 'FaceColor',[(StimulusSet(iStimSet,i)./max(StimulusSet(iStimSet,:))) 0.5 0.5])
                xlabel('Voxel ID','FontSize', fFontSize), ylabel('Response Level','FontSize', fFontSize);
                title ([StimulusSetNames{iStimSet},'Response'], 'FontSize', fFontSize);
                hold on
            end
            if iStimSet == 4
                uicontrol('Style', 'text',...
                    'String', TonotopicNames{iTonoMap},...
                    'FontSize', 25,...
                    'Units','normalized',...
                    'backgroundcolor',[1 1 1],...
                    'ForegroundColor',[0 0 0],...
                    'Position', [0 0.5 1 0.05]);
            end
        end
    end
end

if fESTpCFSel == 1;
    %% Plot all pCFs - all voxels and frequencies
    % Create figure plotting estimated pCF for each magnification domain for
    % each stimulus set with axis limits set by Auditory system min and maximum
    % frequency
    for iTonoMap = 1:length(TonotopicMagnification)
        fTonotopic(iTonoMap) = figure('Color', [1 1 1]);
        set(fTonotopic(iTonoMap),'OuterPosition',FigSize)
        plot (VoxelDistance,TonotopicMap(iTonoMap).pCF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2);
        axis tight, title ([TonotopicNames{iTonoMap},' - Range: ', num2str(LowFreq),' - ',num2str(HighFreq),' kHz'], 'FontSize', 20);
        ylim ([LowFreq HighFreq]);
        hold on
        %     plot ((nNeurons./NeuronsPerVoxel),TonotopicMap(iTonoMap).CF, 'color','b', 'LineWidth', 2);
        for iStimSet = 1:length(StimulusSetNames)
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            %         EstpCF = [zeros(1,xmin), EstpCF(xmin+1:xmax), zeros(1,length(EstpCF) - xmax)];
            scatter (VoxelDistance,EstpCF, 50, SpacingColour{iStimSet}, 'fill');
        end
        set(gca,'YScale','log');
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        fLeg = legend('Actual',...
            [StimulusSetNames{1},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSE(1))],...
            [StimulusSetNames{2},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSE(2))],...
            [StimulusSetNames{3},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSE(3))],...
            [StimulusSetNames{4},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSE(4))],...
            'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
    
    %% Save RMSE to excel file
    T = array2table(vertcat(TonotopicMap.RMSE),'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'},...
        'RowNames',{'FreqDis' 'ERB' 'Log'});
    writetable(T,'RMSE.csv');
end

if fStimulusvsTonotopicSel == 1
    for iStimSet = 1:length(StimulusSetSpacing)
    fStimTono = figure('Color', [1 1 1]);
    set(fStimTono,'OuterPosition',FigSize)
    % create figure to plot stimulus properties
    
        for iTonoMap = 1:length(TonotopicMagnification)
            plot (VoxelDistance,TonotopicMap(iTonoMap).ESTpCF(iStimSet,:));
            hold on
        end
    end
end


end
%% Local Functions
function properties = lcfCreateTonotopicProperties(Fn,TWn,f,nNeurons,G,Gn,Gi,C)
TWn = 0.1; % tuning width at Fn
Fn = 2; % Frequency to set D (point image) to zero (kHz)
% StepPercentage = 0.001;

A = (TWn * Fn) / G(Fn);
B = -(A * (Gn(Fn)));
D = A * (Gn(f)) + B;

f0 = 0.02;
f1 = 20;
x1 = A * (Gn(f0)) + B;
x2 = A * (Gn(f1)) + B;
% so1 = 10^(floor(log10(abs(x1))));
% x1 = so1 * round(x1/so1);
% so2 = 10^(floor(log10(abs(x2))));
% x2 = so2 * round(x2/so2);
x1 = round(x1);
x2 = round(x2);
step = (x2-x1)/nNeurons;
properties.D = x1:step:x2;
properties.CF = Gi(((properties.D) - B)./A);

properties.df = A * G(properties.CF);

properties.p = C*4*properties.CF ./ properties.df; % Filter coefficient of each Voxel % VTuningWidth - population tuning width / populationERB
end

function properties = NOTUSEDlcfCreateTonotopicProperties(LowFreq,HighFreq,nNeurons,C,TonotopicMagnification,TonoMag2Freq,TuningWidth,nr)
%(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
TonotopicFrequencySpacing = linspace(TonotopicMagnification(LowFreq),TonotopicMagnification(HighFreq),nNeurons); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
% Characterisic Frequencies
properties.CF = TonoMag2Freq(TonotopicFrequencySpacing);   % neuron Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
% Tuning Width
% df = TuningWidth(properties.CF); % Neuron frequency tuning width for each Voxel Charactersic Frequency (for each neuron)(Q)

% step = (TonotopicMagnification(HighFreq) - TonotopicMagnification(LowFreq)) /nNeurons;
step = (HighFreq - LowFreq) /nNeurons;
df = diff(properties.CF)./step;
% df = (diff(properties.CF));
% dfnomlise = 1./(diff(TonotopicMagnification(properties.CF))./step);
properties.df = [df(1) df(1:end)];
% dfnomlise = 1./(TonotopicMagnification(max(properties.CF))-TonotopicMagnification(min(properties.CF)));
% dfnormalise1over = 1./(diff(properties.CF)./step);
% dfnormalise1over = [dfnormalise1over(1) dfnormalise1over(1:end)];
% dfnormalise = diff(properties.CF)./step;
% dfnomlise = 1./(diff(TonotopicMagnification(properties.CF))./step);
% dfnormalise = [dfnormalise(1) dfnormalise(1:end)];
% dfnomlise = sum(dfnomlise);
% properties.df = df.*dfnomlise;
% Q = Q(1:end-1)+step./2;

% properties.df = properties.df ./ nr;    % normalise Tuning width at 1 kHz by ERB

properties.p = C*4*properties.CF ./ properties.df; % Filter coefficient of each Voxel % VTuningWidth - population tuning width / populationERB
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
% noise = randn (1,length(CF)).*NoiseLev;
Response = 10*log10(max(Int,1)); % Threshold intensity response and convert to dB
% Response = Response + noise; %% set negative values to 0
% Response = Int; % no threshold
% Response = max(Int,1);
end
function ESTpCF = lcfCalculatepCFEstimation (StimulusFreqs, Response,StimSpacingDomain,StimSpacingDomain2Freq)
%% Estimate population characterstic frequency using weighted average response to each stimulus
% Predefine variable for estimated pCF's
ESTpCF = zeros(length(Response),1);
StimFreqStimDomain = StimSpacingDomain(StimulusFreqs);  % convert stimulus frequencies to stimulus ID - actually stimulus spacing domain
% StimFreqStimDomain = StimulusFreqs ./ max(StimulusFreqs);
for i = 1:length(Response)
    ESTpCF (i) = sum((StimFreqStimDomain'.*Response(:,i)))./sum(Response(:,i)); % need to convert to linear spacing - convert to mag domain then convert back for characteristic frequencies
end
ESTpCF(isnan(ESTpCF)) = 0 ; % remove nans - occurs when voxel doesn't response to a stimulus
% ESTpCF = ESTpCF .* max(StimulusFreqs);
ESTpCF = StimSpacingDomain2Freq(ESTpCF);
%% need to convert back to frequency domain
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
end
function f = lcfInvNDLF10 (nDLF)
%% interpolates to estimate frequency as a function of number of DLFs
fs = 0:0.01:20;
fs = fs.*1000;  % convert from kHz to Hz
a = 0.026;
b = -0.533;
DLF = -(2.^(1-b-a .* sqrt(fs)) .* 5.^(-b-a .* sqrt(fs)) .* (1+a .* sqrt(fs) .* log(10)))./(a.^2 .*(log(2)+log(5)).^2);
f = interp1(DLF,fs,nDLF,'spline');
f = f./1000; % convert from Hz to kHz
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
