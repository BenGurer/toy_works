function data = ModelTonotopticMagnification
% Model of Frequency magnification imaged using fMRI
% effect on stimulus spacing on recorded response compared to underlying
% spacing of neuron/voxel characterisic frequency spacing

% New TO DO
% Use http://www.wolframalpha.com/ to get integral
% Use look up table to get inverse of this – sample for lots of frequencies
% 	Plot as points to identify undersampling
% Use interp1 to interpolate cortical distance value for each Fc
%
% Just focus on ERB and Frequency discrimination cortical spacing
%
% Add a constant ratio to tuning width
% 	Normalise this ratio to the ERB at 1kHz
%
% Change scale levels
% 	Create neurons with pCF - how many neurons per mm in the cortex?
% 	Sample based on constant distance to create voxels – take mean frequency
% 	Implement point spread due to hemodynamic spread


%% TO DO
% Add spread of HMDR - apply to voxel response
% Add noise
% create bandpass noise stimuli
% each voxel contain many tuning curves - loop current code for each voxel
% and take weighted average - then convovle with gaussian Hymodynmaic
% spread response
% 30 random numbers - create coiefient - what is this and why did I say it?
%
% convert SPL to SL for Frequency Discrimination equation
%
% calculate error
% squared error
% (actual2 - estimated2) of each voxel
% mean of voxels
% square root (mean of voxels)

% root mean square error
% for iVMD = 1:length(iVMD)
% for iSMD = 1:length(iSMD)
% for i = 1:nVoxel
% SE(i)=(actual - estimated).^2;
% end
% MSE(iSMD) = mean(SE)
% end
% RMSE(iVMD) = sqrt(MSE(iSMD))
% end

% has frequency discrimation been tested at high frequencies or just
% exprolated

% set xlim of pCF zoom by finding voxels whos pCF frequency matchs
% stimulus frequency limits

% convert voxels to milimeters
% save number of voxels within stimulus set frequency range
% calculate error within stimulus set frequency range

clc; close all

%% Stimuli Magnification domain functions
% magDomain,magScale2F,magFilterWidth
StimMagDomain = {str2func('@LinMagDom'),str2func('@lcfNDLF'),str2func('@lcfNErb'),str2func('@LogMagDom')};
StimMagDomain2Freq = {str2func('@LinMagDom2Freq'),str2func('@lcfInvNDLF'),str2func('@lcfInvNErb'),str2func('@LogMagDom2Freq')};

StimSpacingNames = {'Linear','Frequency Discrimination','Equivalent Rectangular Bandwidth','Logarithmic'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables
% StimSpacingNamesv2 = {'Linear','Logarithmic','ERB','Freq Dis'}; % shorter names if needed due to space

%% Stimuli Properties
StimHighFreq = 8;  % highest frequency to test
StimLowFreq = 0.25; % lowest frequency to test
StimNum = 10;   % number of stimulus
StimdB = 70;
%% Voxel Population Properties
% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.02; % Lowest frequency of the system
Gradient = 30; % length of gradient in mm
nNeurons = 10000; % number of neurons per mm

%% Neuronal population Magnification domain fucntion
% VoxelMagSpace = {str2func('@lcfNDLF'),str2func('@lcfNErb')};
VoxelMagSpace =StimMagDomain;
% VoxelMagSpace2Freq = {str2func('@lcfInvNDLF'),str2func('@lcfInvNErb')};
VoxelMagSpace2Freq = StimMagDomain2Freq;
VoxelFilterWidth = {str2func('@lcfErb'), str2func('@lcfDLF'),str2func('@lcfErb'),str2func('@lcfErb')};
nr = TWnormalise (1,@lcfErb,@lcfDLF); %(f,TWrFun,TWnFun)
% VoxelFilterWidth = {str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb'),str2func('@FreqDisTuningWidth')};
% Spacing Properites
C = 1; % constant scaling for voxel tuning width
NoiseLev = 10; % level of noise in system

%% Imaging Properties
Resolution = 1;  % resolution of imaging in mm
nVoxels = Gradient/Resolution; % number of voxels
NeuronDensity = nNeurons .* Resolution; % number of neurons in each voxel
HDspread = 3; % spread of hemodynamic response in mm

%% Figure Properites
% create colours to use for plot each stimulus/voxel property
SpacingColour = {[1 0 0],[0 0.8 0],[0 0 1],[0 0.8 0.8]};
% use code below to set size of each figure created
% set(f,'OuterPosition',FIG)
FigSize = lcfSetFigSize;
fFontSize = 20;

%% Set up varibles to save too
VMD = cell(length(VoxelMagSpace),1); % Voxel properties - use length for looping
SMD = struct('VoxRes', ones(StimNum,nVoxels), ... % Stimulus propertie - use length for looping
    'VCharaFreq', ones(nVoxels), ...
    'fMRIdata', ones(nVoxels));
Res = cell(length(StimMagDomain),1); % save voxel responses too

%% Create figure to display stimulus properties
fStim = figure('Color', [1 1 1]);
set(fStim,'OuterPosition',FigSize)
% fOneResponse = figure('Color', [1 1 1]);
% set(fOneResponse,'OuterPosition',FIG)
% fOneVoxelResponse = figure('Color', [1 1 1]);
% set(fOneVoxelResponse,'OuterPosition',FIG)
% fOneStim = figure('Color', [1 1 1]);
% set(fOneStim,'OuterPosition',FIG)
% fOnepCF = figure('Color', [1 1 1]);
% set(fOnepCF,'OuterPosition',FIG)
% fOneVoxelHist = figure('Color', [1 1 1]);
% set(fOneVoxelHist,'OuterPosition',FIG)
for iVMD = 1:length(VoxelMagSpace)
    
    %% pCF CREATION
    % Create population charactersic frequency (pCF) Spacing Properites
    data.VoxelProp(iVMD)  = VoxelProperties(LowFreq,HighFreq,nVoxels,C,VoxelMagSpace{iVMD},VoxelMagSpace2Freq{iVMD},VoxelFilterWidth{iVMD},nr);    %(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
    
    %% Voxel Response Loop
    % Create figure to plot stimulus responses
    %     fResponse = figure('Color', [1 1 1]);
    %     set(fResponse,'OuterPosition',FIG)
    %     fVoxelResponse = figure('Color', [1 1 1]);
    %     set(fVoxelResponse,'OuterPosition',FIG)
    %
    for iSMD = 1:length(StimMagDomain)
        
        %% Stimulus Set Creation
        % create stimuli frequencies for each spacing
        StimFreq = StimMagDomain2Freq{iSMD}(linspace(StimMagDomain{iSMD}(StimLowFreq),StimMagDomain{iSMD}(StimHighFreq),StimNum)); % frequencies of stimuli - 3D matrix for band limited noise
        % create stimuli gain structure
        StimLev = ones(size(StimFreq)).*StimdB; % can use to model filter response
        
        %% pCF Responses
        % loop for each stimulus frequency
        for i = 1:length(StimFreq)
            %             % plot to response figure
            %             figure (fResponse)
            % Response from all voxels to stimulus presented- looped for all stimulus properties and looped for all Voxel properties
            SMD.VoxRes(i,:) = VoxelResponse (StimFreq(i),StimLev,data.VoxelProp(iVMD).VFreqC,data.VoxelProp(iVMD).p,NoiseLev);
            SMD.VoxRes(SMD.VoxRes<0) = 0; % remove negative values
            %             subplot(length(VoxelMagSpace)./2,length(StimMagDomain)./2,iSMD); % create and index subplot
            %             area(data.VoxelProp(iVMD).VFreqC,SMD.VoxRes(i,:), 'FaceColor',[(StimFreq(i)./max(StimFreq)) 0.5 0.5]); % area plot of each population response to each stimulus over layed for each stimulus set
            %             % set plot properties
            %             set(gca,'XScale','log');
            %             xlabel('pCF (kHz)', 'FontSize', fFontSize), ylabel('Response Level', 'FontSize', fFontSize);
            %             axis tight, ylim ([0 StimLev(i)+20]), xlim ([LowFreq 20]);
            %             title (StimSpacingNames{iSMD}, 'FontSize', fFontSize);
            %             hold on
            %             % plot stimulus frequencies on graphs for reference
            %             line(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 5)
            %             % add title to centre of figure
            %             if iSMD == 4
            %                 uicontrol('Style', 'text',...
            %                     'String', StimSpacingNames{iVMD},...
            %                     'FontSize', 25,...
            %                     'Units','normalized',...
            %                     'backgroundcolor',[1 1 1],...
            %                     'ForegroundColor',[0 0 0],...
            %                     'Position', [0 0.5 1 0.05]);
            %             end
            
            %             if iSMD == 3 && i == 5 && iVMD == 3
            %                 figure(fOneResponse)
            %                 area(data.VoxelProp(iVMD).VFreqC,SMD.VoxRes(i,:), 'FaceColor',[(StimFreq(i)./max(StimFreq)) 0.5 0.5]); % area plot of each population response to each stimulus over layed for each stimulus set
            %                 % set plot properties
            %
            %                 xlabel('pCF (kHz)', 'FontSize', fFontSize), ylabel('Response Level', 'FontSize', fFontSize);
            %                 axis tight, ylim ([0 StimLev(i)+20]), xlim ([LowFreq 20]);
            %                 set(gca,'XScale','log');
            %                 title ('Response', 'FontSize', fFontSize);
            %                 hold on
            %                 % plot stimulus frequencies on graphs for reference
            %                 line(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 5)
            %
            %                 figure(fOneVoxelResponse)
            %                 area(SMD.VoxRes(i,:), 'FaceColor',[(StimFreq(i)./max(StimFreq)) 0.5 0.5]); % area plot of each population response to each stimulus over layed for each stimulus set
            %                 xlabel('Voxel ID','FontSize', fFontSize), ylabel('Response Level','FontSize', fFontSize);
            %                 title ('Response', 'FontSize', fFontSize);
            %                 hold on
            %
            %                 figure(fOneStim)
            %                 %                 line(repmat(StimFreq(1),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 5)
            %                 line(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 5)
            %                 % set plot properties
            %
            %                 xlabel('pCF (kHz)', 'FontSize', fFontSize), ylabel('Response Level', 'FontSize', fFontSize);
            %                 ylim ([0 StimLev(i)+20]), xlim ([LowFreq 20]);
            %                 set(gca,'XScale','log');
            %                 title ('Stimulus', 'FontSize', fFontSize);
            %
            %
            %                 figure (fOnepCF)
            %                 plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2);
            %                 axis tight, title ('pCF', 'FontSize', fFontSize);
            %                 xlim([0 30]);
            %                 ylim([LowFreq HighFreq]);
            %                 %     set(gca,'YScale','log');
            %                 xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('Population Characteristic Frequency pCF (kHz)', 'FontSize', fFontSize);
            %                 warning('off','MATLAB:legend:IgnoringExtraEntries')
            %                 legend((StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast', 'FontSize', fFontSize);
            %
            %                 figure (fOneVoxelHist)
            % %                 bar(SMD.VoxRes(:,iVMD))
            %                  subplot(length(VoxelMagSpace)./2,length(StimMagDomain)./2,iSMD);
            %                 bar(SMD.VoxRes(:,8))
            %                  axis tight, title ('Voxel Response Profile', 'FontSize', fFontSize);
            %                  xlabel('Stimulus ID', 'FontSize', fFontSize), ylabel('Response level', 'FontSize', fFontSize);
            %             end
        end
        
        %         for i = 1:length(StimFreq)
        %             % plot to response figure
        %             figure (fVoxelResponse)
        %             % Response from all voxels to stimulus presented- looped for all stimulus properties and looped for all Voxel properties
        %             %             SMD.VoxRes(i,:) = VoxelResponse (StimFreq(i),StimLev,data.VoxelProp(iVMD).VFreqC,data.VoxelProp(iVMD).p,NoiseLev);
        %             %             SMD.VoxRes(SMD.VoxRes<0) = 0; % remove negative values
        %             subplot(length(VoxelMagSpace)./2,length(StimMagDomain)./2,iSMD); % create and index subplot
        %             area(SMD.VoxRes(i,:), 'FaceColor',[(StimFreq(i)./max(StimFreq)) 0.5 0.5]); % area plot of each population response to each stimulus over layed for each stimulus set
        %             xlabel('Voxel ID', 'FontSize', fFontSize), ylabel('Response Level', 'FontSize', fFontSize);
        %             title (StimSpacingNames{iSMD}, 'FontSize', fFontSize);
        %             hold on
        %             if iSMD == 4
        %                 uicontrol('Style', 'text',...
        %                     'String', StimSpacingNames{iVMD},...
        %                     'FontSize', 25,...
        %                     'Units','normalized',...
        %                     'backgroundcolor',[1 1 1],...
        %                     'ForegroundColor',[0 0 0],...
        %                     'Position', [0 0.5 1 0.05]);
        %             end
        %             %         if iSMD == 3 && i == 5 && iVMD == 3
        %             %
        %             %         end
        %         end
        %% Stimulus set plot
        % create figure to plot stimulus properties
        if iVMD == 1
            for i = 1:length(StimFreq)
                figure (fStim)% plot to stimulus properties figure
                subplot(length(VoxelMagSpace)./2,length(StimMagDomain)./2,iSMD); % create and index subplot
                plot(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 10); % at each stimulus frequency create lines the level of each stimulus presenation level
                % set plot properties
                set(gca,'XScale','log');
                title (StimSpacingNames{iSMD}, 'FontSize', fFontSize), xlabel('Frequency (kHz)', 'FontSize', fFontSize), ylabel('Presentation Level (dB SPL)', 'FontSize', fFontSize);
                ylim ([0 StimLev(i)]), xlim ([LowFreq 20]);
                hold on
            end
        end
        
        %% pCF ESTIMATION
        % weighted average response for each voxel in stimuli domain (stimuli domain because stimuli ID used for calulation)
        SMD.VCharaFreq = VoxelCharacterFrequency (StimFreq, SMD.VoxRes,StimMagDomain{iSMD},StimMagDomain2Freq{iSMD});
        
        %% fMRI response data - FEATURE TO ADD
        % convolve with Hemodynamic spread
        %         SMD.fMRIdata(iSMD) = fMRIResponse (SMD.VCharaFreq(iSMD),HDspread);
        
        %% Loop Output
        Res{iSMD}= SMD; % save in cell array to allow looping through each stimulus set/properties
        
    end
    
    %% Modelling Results Data
    data.VMD{iVMD,1} = Res; % save in cell array to allow looping through each population characteric frequency (pCF) magnification domain
end

%% pCF Magnifcation
% create figure comparing all pCF magnifcation domains
fpCF = figure('Color', [1 1 1]);
set(fpCF,'OuterPosition',FigSize)
for iVMD = 1:length(VoxelMagSpace)
    figure (fpCF)
    hold on
    plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2);
    axis tight, title ('pCF', 'FontSize', 20);
    xlim([0 30]);
    ylim([LowFreq HighFreq]);
    %     set(gca,'YScale','log');
    xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
    warning('off','MATLAB:legend:IgnoringExtraEntries')
    legend((StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast', 'FontSize', fFontSize);
end
% fpCFfunct = figure('Color', [1 1 1]);
% set(fpCFfunct,'OuterPosition',FIG)
% for iVMD = 1:length(VoxelMagSpace)
%     figure (fpCFfunct)
%     hold on
%     plot (VoxelMagSpace{iVMD}(LowFreq:0.1:HighFreq),(LowFreq:0.1:HighFreq), 'color',SpacingColour{iVMD}, 'LineWidth', 2);
%     axis tight, title ('pCF', 'FontSize', 20);
%     xlim([0 30]);
%     ylim([LowFreq HighFreq]);
% %     set(gca,'YScale','log');
%     xlabel('Frequency Gradient (mm)'), ylabel('Population Characteristic Frequency pCF (kHz)');
%     warning('off','MATLAB:legend:IgnoringExtraEntries')
%     legend((StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
% end


%% pCF Magnifcation Estimation RESULTS - figures to present modelling results

% create figure plotting estimated pCF for each magnification domain for each stimulus set
fTonotopic = figure('Color', [1 1 1]);
set(fTonotopic,'OuterPosition',FigSize)
for iVMD = 1:length(VoxelMagSpace)
    figure (fTonotopic)
    subplot(length(VoxelMagSpace)./2, length(VoxelMagSpace)./2, iVMD);
    plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2);
    axis tight, title (StimSpacingNames{iVMD}, 'FontSize', 20);
    hold on
    for iSMD = 1:length(StimMagDomain)
        scatter (1:nVoxels,data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq, 50, SpacingColour{iSMD}, 'fill')
        set(gca,'YScale','log');
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        warning('off','MATLAB:legend:IgnoringExtraEntries')
        fLeg = legend('pCF',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
end

% % create figure plotting estimated pCF for each magnification domain for each stimulus set with axis limits set by stimulus set min and maximum
% fTonotopicZoom = figure('Color', [1 1 1]);
% set(fTonotopicZoom,'OuterPosition',FIG)
% for iVMD = 1:length(VoxelMagSpace)
%     figure (fTonotopicZoom)
%     subplot(length(VoxelMagSpace)./2, length(VoxelMagSpace)./2, iVMD);
%     plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2);
%     axis tight, title (StimSpacingNames{iVMD}, 'FontSize', 20);
%     xmin = find(data.VoxelProp(iVMD).VFreqC <= StimLowFreq,1);
%     xmax = find(data.VoxelProp(iVMD).VFreqC >= StimHighFreq,1);  % highest frequency to test
%     ylim ([0.2 8]), xlim ([xmin xmax]);
%     hold on
%     for iSMD = 1:length(StimMagDomain)
%         scatter (1:nVoxels,data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq, 50, SpacingColour{iSMD}, 'fill');
%         set(gca,'YScale','log');
%         xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
%         fLeg = legend('Actual',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
%         set(fLeg, 'FontSize', fFontSize./2);
%     end
% end

% create figure plotting estimated pCF for each magnification domain for each stimulus set with axis limits set by stimulus set min and maximum
fTonotopicZoomPrune = figure('Color', [1 1 1]);
set(fTonotopicZoomPrune,'OuterPosition',FigSize)
data.VoxelsWithinStimLim = zeros(1,iVMD);
for iVMD = 1:length(VoxelMagSpace)
    figure (fTonotopicZoomPrune)
    subplot(length(VoxelMagSpace)./2, length(VoxelMagSpace)./2, iVMD);
    plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2);
    axis tight, title (StimSpacingNames{iVMD}, 'FontSize', 20);
    xmin = find(data.VoxelProp(iVMD).VFreqC >= StimLowFreq,1);
    xmax = find(data.VoxelProp(iVMD).VFreqC >= StimHighFreq,1);  % highest frequency to test
    data.VoxelsWithinStimLim (iVMD) = xmax -xmin;
    ylim ([0.2 8]), xlim ([xmin xmax]);
    hold on
    for iSMD = 1:length(StimMagDomain)
        EstpCF = data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq;
        EstpCF = [zeros(1,xmin), EstpCF(xmin+1:xmax)', zeros(1,length(EstpCF) - xmax)];
        scatter (1:nVoxels,EstpCF, 50, SpacingColour{iSMD}, 'fill');
        set(gca,'YScale','log');
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        fLeg = legend('Actual',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
end


data.RMSE = zeros(1,(iVMD));
for iVMD = 1:iVMD
    MSE = zeros(1,(iSMD));
    for iSMD = 1:iSMD
        SE = zeros(iSMD,nVoxels);
        for i = 1:nVoxels
            SE(iSMD,i)=(data.VoxelProp(iVMD).VFreqC(i) - data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq(i)).^2;
        end
        MSE(iSMD) = mean(SE(iSMD,:));
        data.RMSE(iVMD,iSMD) = sqrt(MSE(iSMD));
    end
end
T = array2table(data.RMSE,'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'});
writetable(T,'RMSE.csv');

data.RMSEStimLim = zeros(1,(iVMD));
for iVMD = 1:iVMD
    MSE = zeros(1,(iSMD));
    xmin = find(data.VoxelProp(iVMD).VFreqC <= StimLowFreq,1);
    xmax = find(data.VoxelProp(iVMD).VFreqC >= StimHighFreq,1);  % highest frequency to test
    pCF = data.VoxelProp(iVMD).VFreqC;
    pCF = pCF(xmin+1:xmax);
    data.VoxelsWithinStimLim (iVMD) = xmax -xmin;
    for iSMD = 1:iSMD
        EstpCF = data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq;
        EstpCF = EstpCF(xmin+1:xmax);
        SE = zeros(iSMD,length(pCF));
        for i = 1:length(pCF)
            SE(iSMD,i)= (pCF(i) - EstpCF(i)).^2;
        end
        MSE(iSMD) = mean(SE(iSMD,:));
        data.RMSEStimLim(iVMD,iSMD) = sqrt(MSE(iSMD));
    end
end
T = array2table(data.RMSEStimLim,'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'});
writetable(T,'RMSEStimLim.csv');
end



%% Location funcations to call
function Response = VoxelResponse (StimF,StimLev,VoxelFreqC,p,NoiseLev)
StimInt = 10.^(StimLev/10);
Int = zeros(1,length(VoxelFreqC)); % Response intensity for each Voxel
% loop through all voxels
for I = 1:length(VoxelFreqC)
    for II = 1:length(StimF)
        g = abs((StimF(II)-VoxelFreqC(I))/VoxelFreqC(I));
        fw=(1+(p(I)*g)).*exp(-p(I)*g); %Two parameter roex
        Int(I) = Int(I)+fw*StimInt(II);
    end
end
% NoiseLev = 10;
noise = randn (1,length(VoxelFreqC)).*NoiseLev;
Response = 10*log10(max(Int,1)); % Threshold intensity response and convert to dB
Response = Response + noise; %% set negative values to 0
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
function VoxelProp = VoxelProperties(LowFreq,HighFreq,nVoxels,C,VoxelMagSpace,magScale2F,magFilterWidth,nr)
VoxelProp.VmagScale = linspace(VoxelMagSpace(LowFreq),VoxelMagSpace(HighFreq),nVoxels); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
VoxelProp.VFreqC = magScale2F(VoxelProp.VmagScale);   % Voxel Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
VoxelProp.VTuningWidth = magFilterWidth(VoxelProp.VFreqC); % Voxel frequency tuning width for each Voxel Charactersic Frequency (for each voxel)(Bandwidth in?)
VoxelProp.VTuningWidth = VoxelProp.VTuningWidth ./ nr;
VoxelProp.p = C*4*VoxelProp.VFreqC./VoxelProp.VTuningWidth; % Filter coefficient of each Voxel % VTuningWidth - population tuning width / populationERB
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
function DLF = lcfDLF (f)
% log DLM =(a*sqrt(f)+B)) Weir 1977
% DLM = 10^(a*sqrt(f)+B))
% for 40 dB SL
a = 0.026;
b = -0.533;
DLF = 10.^(a*sqrt(f)+b);

end
function nDLF = lcfNDLF (f)
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
nDLF = -(2.^(1-b-a .* sqrt(f)) .* 5.^(-b-a .* sqrt(f)) .* (1+a .* sqrt(f) .* log(10)))./(a.^2 .* (log(10)).^2);
end

function f = lcfInvNDLF (nDLF)
fs = 0:0.01:20;
fs = fs.*1000;  % convert from kHz to Hz
a = 0.026;
b = -0.533;
DLF = -(2.^(1-b-a .* sqrt(fs)) .* 5.^(-b-a .* sqrt(fs)) .* (1+a .* sqrt(fs) .* log(10)))./(a.^2 .* (log(10)).^2);
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