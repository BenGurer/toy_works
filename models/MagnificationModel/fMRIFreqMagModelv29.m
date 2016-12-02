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
% assign colour for each scaling

clc; close all

%% Neuronal population Magnification domain fucntion
% VoxelMagSpace = {str2func('@lcfNErb'),str2func('@FreqDisMagDom')};
% VoxelMagSpace2Freq = {str2func('@lcfInvNErb'),str2func('@FreqDisMagDom2Fre')};
% VoxelFilterWidth = {str2func('@lcfErb'),str2func('@FreqDisMagDom2Freq')};

%% Get monitor size to make figures full screen
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

FIG = [scnsize(1),...
    scnsize(2),...
    scnsize(3) - edge,...
    scnsize(4)];

% use code below to set size of each figure created
% set(f,'OuterPosition',FIG)



%% Stimuli Magnification domain functions
% magDomain,magScale2F,magFilterWidth
StimMagDomain = {str2func('@LinMagDom'),str2func('@LogMagDom'),str2func('@lcfNErb'),str2func('@FreqDisMagDom')};
StimMagDomain2Freq = {str2func('@LinMagDom2Freq'),str2func('@LogMagDom2Freq'),str2func('@lcfInvNErb'),str2func('@FreqDisMagDom2Freq')};
StimSpacingNames = {'Linear','Logarithmic','Equivalent Rectangular Bandwidth','Frequency Discrimination'};   % names for graph titles for each stimulus spacing - must be in same order as StimMagDomain and StimMagDomain2Freq variables
StimSpacingNamesv2 = {'Linear','Logarithmic','ERB','Freq Dis'};

%% Stimuli Properties
StimHighFreq = 8;  % highest frequency to test
StimLowFreq = 0.25; % lowest frequency to test
StimNum = 10;   % number of stimulus
StimdB = 70;
%% Voxel Population Properties
% Auditory system properties
HighFreq = 20;  % Highest frequency of the system
LowFreq = 0.05; % Lowest frequency of the system
Gradient = 30; % length of gradient in mm
nNeurons = 100000; % number of neurons per mm

%% Neuronal population Magnification domain fucntion
% VoxelMagSpace = {str2func('@lcfNErb')};
% VoxelMagSpace2Freq = {str2func('@lcfInvNErb')};
VoxelFilterWidth = {str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb'),str2func('@lcfErb')};
VoxelMagSpace = StimMagDomain;
VoxelMagSpace2Freq = StimMagDomain2Freq;
% VoxelFilterWidth = StimMagDomain;
% Imaging Properties
Resolution = 1;  % resolution of imaging in mm
nVoxels = Gradient/Resolution; % number of voxels
NeuronDensity = nNeurons .* Resolution; % number of neurons in each voxel
HDspread = 3; % spread of hemodynamic response in mm
% Spacing Properites
C = 1; % constant scaling for voxel tuning width

%% Figure Properites
SpacingColour = {[1 0 0],[0 0.8 0],[0 0 1],[0 0.8 0.8]};


VMD = cell(length(VoxelMagSpace),1);
SMD = struct('VoxRes', ones(StimNum,nVoxels), ...
    'VCharaFreq', ones(nVoxels), ...
    'fMRIdata', ones(nVoxels));
Res = cell(length(StimMagDomain),1);
PlotIndex = [0 4 8 12]

fStim = figure('Color', [1 1 1]);
set(fStim,'OuterPosition',FIG)
for iVMD = 1:length(VoxelMagSpace)
    % Spacing Properites
    data.VoxelProp(iVMD)  = VoxelProperties(LowFreq,HighFreq,nVoxels,C,VoxelMagSpace{iVMD},VoxelMagSpace2Freq{iVMD},VoxelFilterWidth{iVMD});    %     (LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
    fResponse = figure('Color', [1 1 1]);
    set(fResponse,'OuterPosition',FIG)
    for iSMD = 1:length(StimMagDomain)
        %% create stimuli frequencies for each spacing
        StimFreq = StimMagDomain2Freq{iSMD}(linspace(StimMagDomain{iSMD}(StimLowFreq),StimMagDomain{iSMD}(StimHighFreq),StimNum)); % frequencies of stimuli - 3D matrix for band limited noise
        
        %% create stimuli gain structure
        StimLev = ones(size(StimFreq)).*StimdB; % can use to model filter response
        
        %% loop for each voxel
        
        for i = 1:length(StimFreq)
            figure (fResponse)
            SMD.VoxRes(i,:) = VoxelResponse (StimFreq(i),StimLev,data.VoxelProp(iVMD).VFreqC,data.VoxelProp(iVMD).p);
            SMD.VoxRes(SMD.VoxRes<0) = 0;
            subplot(length(VoxelMagSpace)./2,length(StimMagDomain)./2,iSMD);
            area(data.VoxelProp(iVMD).VFreqC,SMD.VoxRes(i,:), 'FaceColor',[(StimFreq(i)./max(StimFreq)) 0.5 0.5]), set(gca,'XScale','log'), xlabel('pCF (kHz)'), ylabel('Response Level'), axis tight, ylim ([0 StimLev(i)+20]), xlim ([0 20]);
            title (StimSpacingNames{iSMD}, 'FontSize', 15)
            hold on
            line(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 5)
            %             if iVMD == 1
            %                title (StimSpacingNamesv2{iSMD}, 'FontSize', 15)
            %             end
            if iSMD == 4
                %                ylabel (StimSpacingNamesv2{iVMD}, 'FontSize', 15)
                uicontrol('Style', 'text',...
                    'String', StimSpacingNames{iVMD},...
                    'FontSize', 25,...
                    'Units','normalized',...
                    'backgroundcolor',[1 1 1],...
                    'ForegroundColor',[0 0 0],...
                    'Position', [0 0.5 1 0.05]);
            end
        end
        if iVMD == 1
            for i = 1:length(StimFreq)
                figure (fStim)
                subplot(length(VoxelMagSpace)./2,length(StimMagDomain)./2,iSMD);
                %             plot(data.VoxelProp(iVMD).VFreqC,SMD.VoxRes(i,:), 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5]), set(gca,'XScale','log'), title (StimSpacingNames{iSMD}), xlabel('Actual Voxel Characterisic Frequency (kHz)'), ylabel('Response Level'), axis tight, ylim ([0 StimLev(i)]);
                %             semilogx(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'LineWidth', 2)
                %             set(gca,'XScale','log'), title (StimSpacingNames{iSMD}), xlabel('Actual Voxel Characterisic Frequency (kHz)'), ylabel('Response Level'), axis tight, ylim ([0 StimLev(i)]), xlim ([0 20]);
                
                plot(repmat(StimFreq(i),1,2),[min(ylim) StimLev(i)], 'color',[(StimFreq(i)./max(StimFreq)) 0.5 0.5], 'LineWidth', 10), set(gca,'XScale','log'), title (StimSpacingNames{iSMD}, 'FontSize', 20), xlabel('Frequency (kHz)'), ylabel('Response Level'), ylim ([0 StimLev(i)]), xlim ([0 20]);
                %         text(iSMD./4, iSMD./4, num2str(StimFreq))
                hold on
            end
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
%     for iVMD = 1:length(VoxelMagSpace)
%         figure (2)
%         for iSMD = 1:length(StimMagDomain)
%             subplot(length(VoxelMagSpace),length(StimMagDomain),iSMD + PlotIndex(iVMD));
%             plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'r'), axis tight
%             hold on
%             scatter (1:nVoxels,data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq)
%             set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('C.Freq (kHz)'),title (StimSpacingNames{iSMD}),legend('Actual','Measured','Location','SouthEast')
%         end
%     end
fpCF = figure('Color', [1 1 1]);
set(fpCF,'OuterPosition',FIG)
for iVMD = 1:length(VoxelMagSpace)
    figure (fpCF)
    hold on
    plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2), axis tight, title ('pCF', 'FontSize', 20)
    set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Population Characteristic Frequency pCF (kHz)'),legend((StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
    end

fTonotopic = figure('Color', [1 1 1]);
set(fTonotopic,'OuterPosition',FIG)
for iVMD = 1:length(VoxelMagSpace)
    figure (fTonotopic)
    subplot(length(VoxelMagSpace)./2, length(VoxelMagSpace)./2, iVMD);
    plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2), axis tight, title (StimSpacingNames{iVMD}, 'FontSize', 20)
    hold on
    for iSMD = 1:length(StimMagDomain)
        scatter (1:nVoxels,data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq, 50, SpacingColour{iSMD}, 'fill')
        
        set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Population Characteristic Frequency pCF (kHz)'),legend('pCF',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
        %  xlabel('Frequency Gradient (mm)'), ylabel('Population Characteristic Frequency pCF (kHz)'),legend('Actual',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast')
        
    end
    %     lsline
end
fTonotopicZoom = figure('Color', [1 1 1]);
set(fTonotopicZoom,'OuterPosition',FIG)
for iVMD = 1:length(VoxelMagSpace)
    figure (fTonotopicZoom)
    subplot(length(VoxelMagSpace)./2, length(VoxelMagSpace)./2, iVMD);
    plot (1:nVoxels,data.VoxelProp(iVMD).VFreqC, 'color',SpacingColour{iVMD}, 'LineWidth', 2), axis tight, title (StimSpacingNames{iVMD}, 'FontSize', 20), ylim ([0.2 8])
    hold on
    for iSMD = 1:length(StimMagDomain)
        scatter (1:nVoxels,data.VMD{iVMD,1}{iSMD ,1}.VCharaFreq, 50, SpacingColour{iSMD}, 'fill')
        set(gca,'YScale','log'), xlabel('Frequency Gradient (mm)'), ylabel('Population Characteristic Frequency pCF (kHz)'),legend('Actual',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast');
        %  xlabel('Frequency Gradient (mm)'), ylabel('Population Characteristic Frequency pCF (kHz)'),legend('Actual',(StimSpacingNames{1}),(StimSpacingNames{2}),(StimSpacingNames{3}),(StimSpacingNames{4}),'Location','SouthEast')
        
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
NoiseLev = 10;
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
function VoxelProp = VoxelProperties(LowFreq,HighFreq,nVoxels,C,magDomain,magScale2F,magFilterWidth)
VoxelProp.VmagScale = linspace(magDomain(LowFreq),magDomain(HighFreq),nVoxels); % Linearly space in magnification scale domain (ERB, FS, LOG, LIN)
VoxelProp.VFreqC = magScale2F(VoxelProp.VmagScale);   % Voxel Characterisic Frequency - Convert maginification scale values to frequency values (kHz)
VoxelProp.VTuningWidth = magFilterWidth(VoxelProp.VFreqC); % Voxel frequency tuning width for each Voxel Charactersic Frequency (for each voxel)(Bandwidth in?)
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