function TempTinEEG
clear all;          % Clear workspace
close all;          % Close all figures
clc;                % Clear command Window
cPath = fullfile('/N','matlab','PhD','TempTin');                    % Di
rcxPath = fullfile('/N','matlab','PhD','TempTin','TempTinEEG.rcx'); % RP circuit path
cd(cPath);          % Change directory


%% Get experiment information
%% Use default parameters?
temp = input('Use default parameters? (y/n):', 's');
switch temp
    case 'y'
        data.expData = defExpInfo;
    case 'n'
        data.expData = getExpInfo;
end

%% Load RP2
RP = Circuit_Loader(rcxPath); % Runs Circuit_Loader

%% Create structure to contain results
data.Results = struct('PerceivedLoudness', double.empty, ...
    'Localize', double.empty, ...
    'Spectral', double.empty, ...
    'Loudness', double.empty, ...
    'Likeness', double.empty, ...
    'tinnitusDur', double.empty, ...
    'InduceTime',double.empty);
%% Run experiment
data.Results.InduceTime = 0;
data.Results.tinnitusDur = 0;
if all(bitget(RP.GetStatus,1:3))
    %% Instructions to name EEG recording data
    if skip == 0
        disp('EEG Acquisiton...Prime EEG Recording')
        disp(['Name file: ', data.expData.filename, 'Pre'])
        
        StartInstructions (data.expData)
        %% Collect Resting State EEG Data - Pre Tinnitus
        EEGcollectionTimerEnd (RP, data.expData);
        
        %% Initial induction of tinnitus
        InduceInstructions (data.expData)
        
        TinIn = 0;
        while (TinIn == 0)
            %% Induce tinnitus
            TinntusInduce (RP, data.expData, data.expData.duration);
            data.Results.InduceTime = data.Results.InduceTime + data.expData.duration;
            %% Check if tinnitus has been induced
            TinIn = checkTinStat (data.expData);
        end
        
        %% Instructions to name EEG recording data
        disp('EEG Acquisiton...Prime EEG Recording')
        disp(['Name file: ', data.expData.filename, 'Post'])
        %% Collect EEG Data
        data.Results.tinnitusDur = EEGcollectionUserEnd (RP, data.expData);
    end
    CharaInstructions (data.expData)
    TinIn = 0;
    while (TinIn == 0)
        %% Induce tinnitus
        TinntusInduce (RP, data.expData, data.expData.duration);
        data.Results.InduceTime = data.Results.InduceTime + data.expData.duration;
        %% Check if tinnitus has been induced
        TinIn = checkTinStat (data.expData);
    end
    
    
    %     if skip == 1;
    %         TinnitusTesterInstructions (data.expData)
    %     end
    %% Tinnitus tester
    
    
    % begin tinnitus characterisation
    response = perceivedloudness (RP, data.expData);
    data.Results.PerceivedLoudness = response.PrecievedLoudness;
    data.Results.InduceTime = data.Results.InduceTime + response.InduceTime;
    
    %% LOCALISE CHARACTERISATION
    response = localize (RP, data.expData);
    data.Results.Localize = response.Localize;
    data.Results.InduceTime = data.Results.InduceTime + response.InduceTime;
    
    %% SPECTRAL CHARACTERISATION
    response = spectral (RP, data.expData);
    data.Results.Spectral = response.Spectral;
    data.Results.InduceTime = data.Results.InduceTime + response.InduceTime;
    
    %% LOUDNESS CHARACTERISATION
    LoudnessInstructions (data.expData);
    % check tinnitus is still present
    TinIn = checkTinStat (data.expData);
    if TinIn ==0
        while (TinIn == 0)
            %% Induce tinnitus
            TinntusInduce (RP, data.expData, data.expData.duration);
            data.Results.InduceTime = data.Results.InduceTime + data.expData.repDuration;
            %% Check if tinnitus has been induced
            TinIn = checkTinStat (data.expData);
        end
        LoudnessInstructions (data.expData)
    end
    response = loudness (RP, data.expData, data.Results.Spectral);
    data.Results.Loudness = response.TinRes;
    data.Results.InduceTime = data.Results.InduceTime + response.InduceTime;
    
    %% LIKENESS CHARACTERISATION
    LikenessInstructions (data.expData);
    % check tinnitus is still present
    if TinIn ==0
        while (TinIn == 0)
            %% Induce tinnitus
            TinntusInduce (RP, data.expData, data.expData.duration);
            data.Results.InduceTime = data.Results.InduceTime + data.expData.repDuration;
            %% Check if tinnitus has been induced
            TinIn = checkTinStat (data.expData);
        end
        LikenessInstructions (data.expData);
    end
    response = likeness (RP, data.expData, data.Results.Spectral, data.Results.Loudness);
    data.Results.Likeness = response.TinRes;
    data.Results.InduceTime = data.Results.InduceTime + response.InduceTime;
    
    
    %% Save the structure 'data' to file using name specified by experimenter
    save(data.expData.filename, 'data');
    
    ThankYouScreen (data.expData);
    
end
end
function RP = Circuit_Loader(varargin)
% CIRCUIT_LOADER Loads a *.rcx circuit onto a RP2, returns ActiveX control object
% RP = CIRCUIT_LOADER(connectionType, deviceNumber, circuitPath)
%     User input require if no inputs to function
%     Options for connectionType are 'GB' and 'USB'
%     circuitPath does not require an extension.
%     connectionType defaults to 'GB'
%     deviceNumber defaults to 1
%     Note: code must be modified to work with non-RP2 devices or with *.rco files

if nargin == 3
    
    connectionType = varargin{1};
    deviceNumber = varargin{2};
    circuitPath = varargin{3};
    
elseif nargin == 1
    
    connectionType = 'USB';
    deviceNumber = 1;
    circuitPath = varargin{1};
    
elseif nargin == 0
    
    % path - set this to wherever the examples are stored
    path = 'C:\Users\beng\Documents\MATLAB\C84PRO';
    
    connectionType = input('Enter the type of connection (USB or GB):  ','s');
    connectionType = upper(connectionType);
    
    % Error check for correct connection
    if ~(strcmp(connectionType,'USB')||strcmp(connectionType,'GB'))
        connectionType = 'GB';
        disp('   Device connection = GB');
    end
    
    deviceNumber = input('Enter the device number:  ');
    
    % Error check for correct device number
    if (~isnumeric(deviceNumber) || deviceNumber < 1)
        deviceNumber = 1;
        disp('   Device number = 1');
    end
    
    % Show available circuits
    disp(' ');
    disp(['path: ' path]);
    dir(path)
    
    circuitPath = input('Enter the name of the circuit:  ','s');
    circuitPath = strcat(path,circuitPath);
    
else
    error('Invalid number of arguments.');
end

% Error check circuit file path
if size(strfind(circuitPath,'.rcx')) == 0
    circuitPath = strcat(circuitPath,'.rcx');
end

% Error check for existing file
fileExists=(exist(circuitPath,'file'));
if fileExists==0
    disp('   File doesnt exist'); return;
end

% Load circuit onto device and run
RP = actxcontrol('RPco.x',[5 5 26 26]);

RP.ConnectRP2(connectionType, deviceNumber); % Connects RP2 via USB or GB given the proper device number
RP.Halt; % Stops any processing chains running on RP2
RP.ClearCOF; % Clears all the buffers and circuits on that RP2
disp(['Loading ' circuitPath]);
RP.LoadCOF(circuitPath); % Loads circuit
RP.Run; % Starts circuit

status=double(RP.GetStatus); % Gets the status
if bitget(status,1)==0; % Checks for connection
    disp('Error connecting to RP2'); return;
elseif bitget(status,2)==0; % Checks for errors in loading circuit
    disp('Error loading circuit'); return;
elseif bitget(status,3)==0 % Checks for errors in running circuit
    disp('Error running circuit'); return;
else
    disp('Circuit loaded and running'); return;
end
end
% Experiemental constants
% Changable parameters of stimulus
% data collection options
function [data] = defExpInfo
% create structure to contain the information needed
data = struct('filename', char.empty, ...
    'SF', double.empty, ...
    'Cf', double.empty, ...
    'location', double.empty, ...
    'filter',double.empty, ...
    'BWERB',double.empty, ...
    'duration',double.empty,...
    'repDuration',double.empty,...
    'acqTimeMin',double.empty,...
    'EyesOpenAcqTimeMin',double.empty,...
    'EyesClosedAcqTimeMin',double.empty,...
    'scnsize',double.empty,...
    'edge',double.empty,...
    'stim_dB', double.empty, ...
    'NomLev', double.empty,...
    'HB7Gain',double.empty);

%% Sampling frequency
data.SF = 50000;

%% Set Level
%     data.NomLev = 109.0397; % HD600
data.NomLev = 113; % HD280
data.HB7Gain = -15;
data.stim_dB = 60;

%% Presentation time
%     data.duration = 60;
%     data.repDuration = 10;

data.duration = 0.01;
data.repDuration = 0.01;

%% EEG acquisition time
%     data.EyesOpenAcqTimeMin = 5;
%     data.EyesClosedAcqTimeMin = 5;
data.EEGAcqTimeMin  = 0.01;

%% Set figure size

f = figure ('Visible','off',...
    'Menu', 'None', ...
    'Color', [0.5, 0.5, 0.5]);

set(0,'Units','pixels') ;
p = get(0, 'MonitorPositions');
H = size (p);

% switch H(1)
%     case 1
data.scnsize = p(1,:);
%     case 2
%         m1 = p(1,:);
%         m2 = p(2,:);
%         x = m2(3) - m1(3);
%         data.scnsize = [m2(1), m2(2) - (m2(4).*0.165), x, m2(4)];
% end

position = get(f,'Position');
outerpos = get(f,'OuterPosition');
borders = outerpos - position;
data.edge = -borders(1)/2;

close (f)
%% Request the name of the filename to save results and experiment data to
data.filename = input('Enter the filename: ', 's');

% Check if information has been entered. Ask until data given.
check = 0;
while (check == 0)
    done = 0;
    while( done==0 )
        if  ~isempty(data.filename)
            valuefilenameCheck = 1;
        else
            disp ('Please enter filename')
            valuefilenameCheck = 0;
        end
        
        if( valuefilenameCheck == 1 )
            done = 1;
        else
            data.filename = input('Enter the filename: ', 's');
        end
    end
    % Check if the filename contains spaces. Ask until there isn't.
    done = 0;
    while( done==0 )
        if sum(isspace(data.filename)) == 0;
            filenameCheck = 1;
        else
            disp ('Filename cannot contain spaces')
            filenameCheck = 0;
        end
        
        if( filenameCheck == 1 )
            done = 1;
        else
            data.filename = input('Enter the filename: ', 's');
        end
    end
    check = 1;
end
%% Cut-off frequency of inducing stimulus in kHz

data.Cf = 1;

%% Presenation ear for inducing stimulus
%         1 = left
%         2 = right
data.location = 1;
%                 data.location = 2;


%% Filter type for inducing stimulus
% 'Lowpass (1), Highpass (2) or Notch (3)?

data.filter  = 1;
%                 data.filter  = 2;
%                 data.filter  = 3;

if (data.filter == 3)
    prompt = input('Notch distance from filter cut-off (ERB):');
    data.BWERB = input(prompt);
else
    data.BWERB = 0;
end

end
function [data] = getExpInfo
% create structure to contain the information needed
data = struct('filename', char.empty, ...
    'SF', double.empty, ...
    'Cf', double.empty, ...
    'location', double.empty, ...
    'filter',double.empty, ...
    'BWERB',double.empty, ...
    'duration',double.empty,...
    'repDuration',double.empty,...
    'acqTimeMin',double.empty,...
    'EyesOpenAcqTimeMin',double.empty,...
    'EyesClosedAcqTimeMin',double.empty,...
    'scnsize',double.empty,...
    'edge',double.empty,...
    'stim_dB', double.empty, ...
    'NomLev', double.empty,...
    'HB7Gain',double.empty);


%% ASK EXPERIMENTER TO INPUT REQUIRED INFORMATION.
% Input through command line
% error checking will break if user does both errors tested for -
% could be fixed by making one while loop to encompass both error
% checks

%% Sampling frequency
data.SF = 50000;

% Set Level
%     data.NomLev = 109.0397; % HD600
data.NomLev = 113; % HD280
data.HB7Gain = -15;
data.stim_dB = 60;

%% Set figure size

f = figure ('Visible','off',...
    'Menu', 'None', ...
    'Color', [0.5, 0.5, 0.5]);

set(0,'Units','pixels');
p = get(0, 'MonitorPositions');
H = size (p);

switch H(1)
    case 1
        data.scnsize = p(1,:);
    case 2
        m1 = p(1,:);
        m2 = p(2,:);
        x = m2(3) - m1(3);
        data.scnsize = [m2(1), m2(2) - (m2(4).*0.165), x, m2(4)];
end

position = get(f,'Position');
outerpos = get(f,'OuterPosition');
borders = outerpos - position;
data.edge = -borders(1)/2;

close (f)
%% Request the name of the filename to save results and experiment data to
data.filename = input('Enter the filename: ', 's');

% Check if information has been entered. Ask until data given.
check = 0;
while (check == 0)
    done = 0;
    while( done==0 )
        if  ~isempty(data.filename)
            valuefilenameCheck = 1;
        else
            disp ('Please enter filename')
            valuefilenameCheck = 0;
        end
        
        if( valuefilenameCheck == 1 )
            done = 1;
        else
            data.filename = input('Enter the filename: ', 's');
        end
    end
    % Check if the filename contains spaces. Ask until there isn't.
    done = 0;
    while( done==0 )
        if sum(isspace(data.filename)) == 0;
            filenameCheck = 1;
        else
            disp ('Filename cannot contain spaces')
            filenameCheck = 0;
        end
        
        if( filenameCheck == 1 )
            done = 1;
        else
            data.filename = input('Enter the filename: ', 's');
        end
    end
    check = 1;
end
%% Request cut-off frequency of inducing stimulus
prompt = 'Enter Cut-off Frequency (kHz):';
data.Cf = input(prompt);

%% Request presenation ear for inducing stimulus
temp = input('Filter left (l) or right (r) ear?:', 's');

switch temp
    case 'l'
        data.location = 1;
    case 'r'
        data.location = 2;
end

%% Request filter type for inducing stimulus
temp = input('Lowpass (l), Highpass (h) or Notch (n)?:', 's');
switch temp
    case 'l'
        data.filter  = 1;
    case 'h'
        data.filter  = 2;
    case 'n'
        data.filter  = 3;
end

if (data.filter == 3)
    temp = input('Notch distance from filter cut-off (ERB):');
    data.BWERB = temp;
else
    data.BWERB = 0;
end

% Presentation time
temp = input('Duration of tinnitus induction (mins):');
data.duration = temp;

% Repeat presentation time
temp = input('Duration of repeat tinnitus induction (mins):');
data.repDuration = temp;


%% EEG acquisition time
temp = input('Duration of EEG acquisition (mins):');
data.EEGAcqTimeMin  = temp;

end
function [data] = genExpPar (n, nr)
%% Function to create psuedorandom order to present phonemes
% input
% n = number of frequencies to test
% nr = number of times to present each frequencies
% output =
%       data.TinTestFreqOrder = trial order
%       data.iTrial = total number of trials

%% Create structure to save number of trials and order of trials too
data = struct('iTrial', double.empty, ...
    'phonemeOrder', double.empty);


%% create matrix with number of frequencies required and save total number
% of trails (iTrial) to use throughout the function
a = [1:n]';
b = repmat(a,nr, 1);
data.iTrial = (length(b));

%% randperm to create random matrix to randomally index
ix = randperm(data.iTrial);

%% index phoneme matrix randomally to create phoneme presentation order
data.TinTestFreqOrder = b(ix);

end
%% Create Tinnitus Inducing Stimuli
function TinStimCreator (RP, expData, duration)
disp('Inducing Tinnitus')
%% expData struc
%     	data = struct('filename', char.empty, ...
%         'SF', double.empty, ...
% 		'Cf', double.empty, ...
%         'location', double.empty, ...
% 		'filter',double.empty, ...
%         'DERB',double.empty, ...
%         'duration',double.empty,...
%         'stim_dB', double.empty, ...
%         'NomLev', double.empty,...
%         'HB7Gain',double.empty);
%
durationMins = (duration .* 60000) ./3;
MinWarning =  ((duration - 1) .* 60000) ./3;

BuffDur = durationMins + 100;

ERBcomp = lcfLee(expData.Cf,expData.SF);
stim_amp=10^((expData.stim_dB-expData.NomLev-ERBcomp-expData.HB7Gain)/20);

RP.SetTagVal('TinAmp',stim_amp);

%% Generate sigal
soundsource = generateTIstim (expData.Cf, expData.BWERB, expData.filter, expData.location, 10, expData.SF);

%% Set TDT parameters
RP.SetTagVal('OutputChannel',0);
RP.SetTagVal('SamLength',soundsource.SampLength);
RP.SetTagVal('Duration',durationMins);
RP.SetTagVal('MinWarnTimer',MinWarning);
RP.SetTagVal('BuffDur',BuffDur);
RP.WriteTagV('LeftNoise',0,soundsource.left);
RP.WriteTagV('RightNoise',0,soundsource.right);

% wait to allow for loading
pause (0.1)

disp('Presenting stimuli')

%% Start presentation
RP.SoftTrg(3);

% wait for presentation to finish
active = logical(RP.GetTagVal('active_trig'));
while active
    active = logical(RP.GetTagVal('active_trig'));
    %     warning = logical(RP.GetTagVal('MinWarning'));
    %     if (warning == 0)
    %         disp('One minute warning')
    %     end
    pause(0.01)
end

%% End presentation

close (f)

disp('Presentation complete')

    function Lee = lcfLee(F,SF)
        % Lee = level of ee-noise within 1 ERB around F of wideband ee-noise relative to overall level;
        F1 = lcfInvNErb(lcfNErb(F)-0.5); F2 = lcfInvNErb(lcfNErb(F)+0.5);
        Lee = 10*log10(lcfIntErb(F1,F2)/lcfIntErb(0,SF/2));
    end
    function I = lcfIntErb(F1,F2)
        % ***** lcfIntErb *****
        % This calculates the integral of 1./erb(f) between F1 and F2;
        A = 24.7/1000; B = 4.37;
        I = (log(A*(B*F2+1))-log(A*(B*F1+1)))/(A*B);
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
    function soundsource = generateTIstim (Fc, BWERB, FilterMode, location, BuffLength, sampclock)
        %%
        %% Set parameters
        % sampclock     % sampling rate (Hz)
        % BuffLength     % duration of buffer in seconds (s)
        % Fs            % Cut off frequency of filters (kHz)
        % BWERB          % Bandwidth of notch filter in ERB (ERB)
        
        %% calculate paramitors
        N = (BuffLength * sampclock);       % number of samples (samples)
        ST = 1/(sampclock/1000);            % length of each sample (ms) (clock in kHz)
        DF = 1/(ST*N);                      % frequency of each sample - how often/long expressed as Hz (kHz)
        frq = DF*(1:N/2);                   % Highest frequencies represented (nyquist frequency) (kHz)
        
        %% Create sound sources
        soundsource = lcfMakeFNoise(N, DF, frq, Fc, FilterMode, location, BWERB);
        soundsource.SampLength = N;
        plotsound(N, Fc, soundsource.left, soundsource.right, frq, location)
        
    end
    function data = lcfMakeFNoise(N, DF, frq, Fc, FMode, location, BWERB)
        % ******************** lcfMakeFNoise ********************
        %% Set parameters
        % N =           Length of buffer in samples
        % DF =          frequency of each sample - how often/long expressed as Hz (kHz)
        % frq =         Highest frequencies represented (nyquist frequency) (kHz)
        % Fc =          Cut off frequency of filter (kHz)
        % BWERB =       Bandwidth of notch filter in ERB (ERB)
        % FMode =       Select filter type - 1 = LPF, 2 = HPF, 3 = Notch, 0 = no filter
        %%
        notch = ones(1,N/2);                % create brick wall filter
        if FMode==1
            notch(round(Fc/DF):end) = 0;    % low pass filter
        elseif FMode==2
            notch(1:round(Fc/DF)) = 0;      % high pass filter
        elseif FMode==3
            NErbFs = lcfNErb(Fc);               % convert Hz to ERB
            Fc = lcfInvNErb(NErbFs+(BWERB/2));       % add distance to Centre frequency in ERB then convert to Hz
            Fclow = lcfInvNErb(NErbFs-(BWERB/2));    % low cutoff of notch
            
            notch(round(Fclow/DF):round(Fc/DF)) = 0;  % notch filter (DERB wide around Fc)
        end
        
        lev = -10*log10(lcfErb(frq));       % create equivilent values for each energy in each ERB based on highest frequency possible to represent
        filter = 10.^(lev/20);              % create filter to equally weight energy in each ERB
        
        noise = randn(1,N);                 % create noise length of sample N
        noise = real(ifft([filter fliplr(filter)].*fft(noise))); % filter with equal power per ERB filter
        noise = noise/sqrt(mean(noise.^2)); % normalise mean level of total signal
        fNoise = real(ifft([notch fliplr(notch)].*fft(noise))); % filter noise with cliff edge filter (filter already in frequency domain)
        
        switch location
            case 1
                data.left = fNoise;
                data.right = noise;
            case 2
                data.left = noise;
                data.right = fNoise;
        end
    end
    function data = plotsound(N, Fc, left, right, frq, location)
        %% FFT sound to plot
        leftfspec = fft(left);
        rightfspec = fft(right);
        
        % Plot sound
        % Create figure
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        set(0,'Units','pixels')
        p = get(0, 'MonitorPositions');
        H = size (p);
        
        switch H(1)
            case 1
                scnsize = p(1,:);
            case 2
                %                 scnsize = p(2,:);
                scnsize = p(1,:);
        end
        
        position = get(f,'Position');
        outerpos = get(f,'OuterPosition');
        borders = outerpos - position;
        edge = -borders(1)/2;
        
        pos =   [scnsize(1),...
            scnsize(4)./2,...
            (scnsize(3)./2) - edge,...
            scnsize(4)./2];
        
        set(f,'OuterPosition',pos)
        
        % Create axes
        axes1 = axes('Parent',f,'YGrid','on',...
            'XScale','linear',...
            'XMinorTick','on',...
            'XMinorGrid','on',...
            'XGrid','on');
        xlim(axes1,[0 20]);
        ylim(axes1,[-125 100]);
        hold(axes1,'all');
        switch location
            case 1
                plot(axes1, frq, (20*log10(abs(rightfspec(1:N/2)+0.000001))),'g-', frq, (20*log10(abs(leftfspec(1:N/2)+0.000001))),'b-')
                legend('Right', 'Left')
            case 2
                plot(axes1, frq, (20*log10(abs(leftfspec(1:N/2)+0.000001))),'g-', frq, (20*log10(abs(rightfspec(1:N/2)+0.000001))),'b-')
                legend('Left', 'Right')
        end
        
        line([Fc Fc],ylim,'Color','k','LineStyle','--')
        xlabel('Frequency (kHz)');
        ylabel('Amplitude (RMS) dB SPL');
        
        set (f, 'Visible', 'on');
        shg
    end
end
function TinStimPlayer (RP, expData, duration)
disp('Inducing Tinnitus')
%% expData struc
%     	data = struct('filename', char.empty, ...
%         'SF', double.empty, ...
% 		'Cf', double.empty, ...
%         'location', double.empty, ...
% 		'filter',double.empty, ...
%         'DERB',double.empty, ...
%         'duration',double.empty,...
%         'stim_dB', double.empty, ...
%         'NomLev', double.empty,...
%         'HB7Gain',double.empty);
%
durationMins = (duration .* 60000) ./3;
MinWarning =  ((duration - 1) .* 60000) ./3;

BuffDur = durationMins + 100;

ERBcomp = lcfLee(expData.Cf,expData.SF);
stim_amp=10^((expData.stim_dB-expData.NomLev-ERBcomp-expData.HB7Gain)/20);

RP.SetTagVal('TinAmp',stim_amp);

%% Generate sigal
soundsource = generateTIstim (expData.Cf, expData.BWERB, expData.filter, expData.location, 10, expData.SF);

%% Set TDT parameters
RP.SetTagVal('OutputChannel',0);
RP.SetTagVal('SamLength',soundsource.SampLength);
RP.SetTagVal('Duration',durationMins);
RP.SetTagVal('MinWarnTimer',MinWarning);
RP.SetTagVal('BuffDur',BuffDur);
RP.WriteTagV('LeftNoise',0,soundsource.left);
RP.WriteTagV('RightNoise',0,soundsource.right);

% wait to allow for loading
pause (0.1)

disp('Presenting stimuli')

%% Start presentation
RP.SoftTrg(3);

% wait for presentation to finish
active = logical(RP.GetTagVal('active_trig'));
while active
    active = logical(RP.GetTagVal('active_trig'));
    %     warning = logical(RP.GetTagVal('MinWarning'));
    %     if (warning == 0)
    %         disp('One minute warning')
    %     end
    pause(0.01)
end

%% End presentation

close (f)

disp('Presentation complete')

    function Lee = lcfLee(F,SF)
        % Lee = level of ee-noise within 1 ERB around F of wideband ee-noise relative to overall level;
        F1 = lcfInvNErb(lcfNErb(F)-0.5); F2 = lcfInvNErb(lcfNErb(F)+0.5);
        Lee = 10*log10(lcfIntErb(F1,F2)/lcfIntErb(0,SF/2));
    end
    function I = lcfIntErb(F1,F2)
        % ***** lcfIntErb *****
        % This calculates the integral of 1./erb(f) between F1 and F2;
        A = 24.7/1000; B = 4.37;
        I = (log(A*(B*F2+1))-log(A*(B*F1+1)))/(A*B);
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
    function soundsource = generateTIstim (Fc, BWERB, FilterMode, location, BuffLength, sampclock)
        %%
        %% Set parameters
        % sampclock     % sampling rate (Hz)
        % BuffLength     % duration of buffer in seconds (s)
        % Fs            % Cut off frequency of filters (kHz)
        % BWERB          % Bandwidth of notch filter in ERB (ERB)
        
        %% calculate paramitors
        N = (BuffLength * sampclock);       % number of samples (samples)
        ST = 1/(sampclock/1000);            % length of each sample (ms) (clock in kHz)
        DF = 1/(ST*N);                      % frequency of each sample - how often/long expressed as Hz (kHz)
        frq = DF*(1:N/2);                   % Highest frequencies represented (nyquist frequency) (kHz)
        
        %% Create sound sources
        soundsource = lcfMakeFNoise(N, DF, frq, Fc, FilterMode, location, BWERB);
        soundsource.SampLength = N;
        plotsound(N, Fc, soundsource.left, soundsource.right, frq, location)
        
    end
    function data = lcfMakeFNoise(N, DF, frq, Fc, FMode, location, BWERB)
        % ******************** lcfMakeFNoise ********************
        %% Set parameters
        % N =           Length of buffer in samples
        % DF =          frequency of each sample - how often/long expressed as Hz (kHz)
        % frq =         Highest frequencies represented (nyquist frequency) (kHz)
        % Fc =          Cut off frequency of filter (kHz)
        % BWERB =       Bandwidth of notch filter in ERB (ERB)
        % FMode =       Select filter type - 1 = LPF, 2 = HPF, 3 = Notch, 0 = no filter
        %%
        notch = ones(1,N/2);                % create brick wall filter
        if FMode==1
            notch(round(Fc/DF):end) = 0;    % low pass filter
        elseif FMode==2
            notch(1:round(Fc/DF)) = 0;      % high pass filter
        elseif FMode==3
            NErbFs = lcfNErb(Fc);               % convert Hz to ERB
            Fc = lcfInvNErb(NErbFs+(BWERB/2));       % add distance to Centre frequency in ERB then convert to Hz
            Fclow = lcfInvNErb(NErbFs-(BWERB/2));    % low cutoff of notch
            
            notch(round(Fclow/DF):round(Fc/DF)) = 0;  % notch filter (DERB wide around Fc)
        end
        
        lev = -10*log10(lcfErb(frq));       % create equivilent values for each energy in each ERB based on highest frequency possible to represent
        filter = 10.^(lev/20);              % create filter to equally weight energy in each ERB
        
        noise = randn(1,N);                 % create noise length of sample N
        noise = real(ifft([filter fliplr(filter)].*fft(noise))); % filter with equal power per ERB filter
        noise = noise/sqrt(mean(noise.^2)); % normalise mean level of total signal
        fNoise = real(ifft([notch fliplr(notch)].*fft(noise))); % filter noise with cliff edge filter (filter already in frequency domain)
        
        switch location
            case 1
                data.left = fNoise;
                data.right = noise;
            case 2
                data.left = noise;
                data.right = fNoise;
        end
    end
    function data = plotsound(N, Fc, left, right, frq, location)
        %% FFT sound to plot
        leftfspec = fft(left);
        rightfspec = fft(right);
        
        % Plot sound
        % Create figure
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        set(0,'Units','pixels')
        p = get(0, 'MonitorPositions');
        H = size (p);
        
        switch H(1)
            case 1
                scnsize = p(1,:);
            case 2
                %                 scnsize = p(2,:);
                scnsize = p(1,:);
        end
        
        position = get(f,'Position');
        outerpos = get(f,'OuterPosition');
        borders = outerpos - position;
        edge = -borders(1)/2;
        
        pos =   [scnsize(1),...
            scnsize(4)./2,...
            (scnsize(3)./2) - edge,...
            scnsize(4)./2];
        
        set(f,'OuterPosition',pos)
        
        % Create axes
        axes1 = axes('Parent',f,'YGrid','on',...
            'XScale','linear',...
            'XMinorTick','on',...
            'XMinorGrid','on',...
            'XGrid','on');
        xlim(axes1,[0 20]);
        ylim(axes1,[-125 100]);
        hold(axes1,'all');
        switch location
            case 1
                plot(axes1, frq, (20*log10(abs(rightfspec(1:N/2)+0.000001))),'g-', frq, (20*log10(abs(leftfspec(1:N/2)+0.000001))),'b-')
                legend('Right', 'Left')
            case 2
                plot(axes1, frq, (20*log10(abs(leftfspec(1:N/2)+0.000001))),'g-', frq, (20*log10(abs(rightfspec(1:N/2)+0.000001))),'b-')
                legend('Left', 'Right')
        end
        
        line([Fc Fc],ylim,'Color','k','LineStyle','--')
        xlabel('Frequency (kHz)');
        ylabel('Amplitude (RMS) dB SPL');
        
        set (f, 'Visible', 'on');
        shg
    end
end
%% Check for Tinnitus
function TinStat = checkTinStat (expData)
disp('Tinnitus characterisation ... Check tinnitus status')

%% Load GUI for participant response and parameter updates
[f,df] = gui (expData);
uiwait(f)
TinStat = getappdata(df,'TinStat');

close (df)

    function [f, df] = gui (expData)
        
        %% Data transfer figure
        % Create a figure to save data to
        df = figure ('Visible','off',...
            'Menu', 'None');
        
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGuiLayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        setappdata(df,'TinStat',[]);
        
        %% Create GUI elements
        % Create push button
        YesButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Yes',...
            'Position', pos.YesBTN,...
            'FontSize', 36,...
            'FontUnits', 'normalized',...
            'FontWeight','bold',...
            'Callback', @Yes);
        
        NoButton = uicontrol(f,'Style', 'pushbutton', 'String', 'No',...
            'Position', pos.NoBTN,...
            'FontSize', 36,...
            'FontUnits', 'normalized',...
            'FontWeight','bold',...
            'backgroundcolor',[1 0 0],...
            'ForegroundColor',[1 1 1],...
            'Callback', @No);
        
        % Add a text uicontrol to label the slider.
        txt = uicontrol(f,'Style','text',...
            'Position',pos.TXT,...
            'String','Are you experiencing tinnitus?',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 72,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        saveas(gcf,'checkTinStat.png')
        
        function Yes(src, e)
            
            setappdata(df,'TinStat',1);
            close (f);
        end
        function No(src, e)
            
            setappdata(df,'TinStat',0);
            close (f);
        end
        function keyPress(src, e)
            
            
            switch e.Key
                
                case 'space'
                    Yes(YesButton, []);
                    
                case 'a'
                    No(NoButton, []);
                    
            end
            
        end
        
        function data = genGuiLayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.5;
            txtSizey = scnsize(4).*0.5;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.YesBTN = [(scnsize(3).*0.35) - (butSizeX./2),...
                (scnsize(4).*0.25) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            data.NoBTN = [(scnsize(3).*0.65) - (butSizeX./2),...
                (scnsize(4).*0.25) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            data.TXT = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.55) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
        end
    end
end
%% EEG Acquisition
function EEGcollectionTimerEnd (RP, expData)

EEGAcquisition = (expData.EEGAcqTimeMin .* 60000) ./3;
MinWarning =  ((expData.EEGAcqTimeMin - 1) .* 60000) ./3;
RP.SetTagVal('MinWarnTimer',MinWarning);
RP.SetTagVal('EEGTimer',EEGAcquisition);

%% Prompt
disp('EEG Acquisiton Ready...Start EEG Recording')
prompt = 'Hit Enter to Begin';
str = input(prompt);

disp('EEG Acquisiton')

[f] = EEGcollectionGUI (RP, expData);
% wait for Acquisiton to finish
active = logical(RP.GetTagVal('EEGActive'));
while active
    active = logical(RP.GetTagVal('EEGActive'));
    %     warning = logical(RP.GetTagVal('MinWarning'));
    %     if (warning ==0)
    %         disp('One minute warning')
    %     end
    pause(0.01)
end
close (f)
[f] = EEGcollectionGUIEXIT (RP, expData);
drawnow

disp('EEG Acquisiton Complete...End EEG Recording')
prompt = 'Hit Enter to Continue';
str = input(prompt);

close (f)

    function [f] = EEGcollectionGUI (RP, expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGUILayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        %% Create GUI elements
        % 		Text
        FontSize = 60;
        
        txt1 = uicontrol(f,'Style','text',...
            'Position',pos.TXTCross,...
            'String','+',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', FontSize,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        %     Create invisible select button
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Continue',...
            'Visible', 'off',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% EEG collection
        RP.SoftTrg(4);
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        saveas(gcf,'EEGcollection.png')
        
        function select(src, e)
            close (f);
        end
        
        function keyPress(src, e)
            
            switch e.Key
                case 'space'
                    select(SelectButton, []);
                    
            end
            
        end
        
        function data = genGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.3;
            txtSizey = scnsize(4).*0.3;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.TXTCross = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.5) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
        end
        % end of functions within GUI
    end
    function [f] = EEGcollectionGUIEXIT (RP, expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGUILayout (f, expData);
        % 		Text
        FontSize = 60;
        %% Create text
        txtCloseEyes = uicontrol(f,'Style','text',...
            'Position',pos.TXTinst,...
            'String','Thank you, please wait...',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', FontSize,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        function data = genGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.3;
            txtSizey = scnsize(4).*0.3;
            txtINSTSizex = scnsize(3).*0.7;
            txtINSTSizey = scnsize(4).*0.7;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.TXTCross = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.5) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.TXTinst = [(scnsize(3).*0.5) - (txtINSTSizex./2),...
                (scnsize(4).*0.5) - (txtINSTSizey./2),...
                txtINSTSizex,...
                txtINSTSizey];
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
        end
        
        % end of functions within GUI
    end
end
function data = EEGcollectionUserEnd (RP, expData)

%% output = data = duration of induced tinnitus

% Prompt
disp('EEG Acquisiton Ready...Start EEG Recording')
prompt = 'Hit Enter to Begin';
str = input(prompt);

disp('EEG Acquisiton')
t_start = clock;
[f] = EEGcollectionGUI (RP, expData);
uiwait(f)
t_stop = clock;
data = etime(t_stop, t_start);
fprintf ( 1, 'Tinnitus dutation = %f\n', data);

[f] = EEGcollectionGUIEXIT (RP, expData);
drawnow

disp('EEG Acquisiton Complete...End EEG Recording')
prompt = 'Hit Enter to Continue';
str = input(prompt);

close (f)

    function [f] = EEGcollectionGUI (RP, expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGUILayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        %% Create GUI elements
        % 		Text
        FontSize = 60;
        
        txt1 = uicontrol(f,'Style','text',...
            'Position',pos.TXTCross,...
            'String','+',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', FontSize,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        %     Create invisible select button
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Continue',...
            'Visible', 'off',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% EEG collection
        RP.SoftTrg(4);
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        
        
        function select(src, e)
            close (f);
        end
        
        function keyPress(src, e)
            
            switch e.Key
                case 'space'
                    select(SelectButton, []);
                    
            end
            
        end
        
        function data = genGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.3;
            txtSizey = scnsize(4).*0.3;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.TXTCross = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.5) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
        end
        % end of functions within GUI
    end
    function [f] = EEGcollectionGUIEXIT (RP, expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGUILayout (f, expData);
        % 		Text
        FontSize = 60;
        %% Create text
        txtCloseEyes = uicontrol(f,'Style','text',...
            'Position',pos.TXTinst,...
            'String','Thank you, please wait...',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', FontSize,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        function data = genGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.3;
            txtSizey = scnsize(4).*0.3;
            txtINSTSizex = scnsize(3).*0.7;
            txtINSTSizey = scnsize(4).*0.7;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.TXTCross = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.5) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.TXTinst = [(scnsize(3).*0.5) - (txtINSTSizex./2),...
                (scnsize(4).*0.5) - (txtINSTSizey./2),...
                txtINSTSizex,...
                txtINSTSizey];
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
        end
        
        % end of functions within GUI
    end
end
%% Create GUI's to rate tinnitus
function data = perceivedloudness (RP, expData)
disp('Tinnitus characterisation ... Perceived Loudness')

done = 0;
data.PrecievedLoudness = NaN (1);
data.InduceTime = 0;
while (done == 0)
    
    
    %% Load GUI for participant response and parameter updates
    [f,df] = PerceivedLoudnessGUI (expData);
    
    uiwait(f)
    TinStat = getappdata(df,'TinStat');
    
    if (TinStat == 0)
        %% run inducing stimuli
        disp ('Tinnitus no longer present... Inducing tinnitus')
        q = i;
        TinIn = 0;
        while (TinIn == 0)
            %             duration = 10; % Time in minutes to present tinnitus inducing stimulus
            TinntusInduce (RP, expData, expData.repDuration);
            data.InduceTime = data.InduceTime + expData.repDuration;
            %% Check if tinnitus has been induced
            % Load yes/no GUI
            TinIn = checkTinStat (expData);
            if (TinIn == 1)
                break
            end
        end
    elseif (TinStat == 1)
        disp ('Tinnitus still present... continue characterisation')
        
        data.PrecievedLoudness = getappdata(df,'Response');
        done = 1;
    end
    
end

close (df)

    function [f,df] = PerceivedLoudnessGUI (expData)
        
        %% Data transfer figure
        % Create a figure to save data to
        df = figure ('Visible','off',...
            'Menu', 'None');
        
        set(0,'Units','pixels')
        scnsize = get(0,'ScreenSize');
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGuiLayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        setappdata(df,'TinStat',1);
        
        %% Create GUI elements
        % Create slider
        sld = uicontrol(f,'Style', 'slider',...
            'Min',0,'Max',100,'Value',0,...
            'backgroundcolor',[0.7 0.7 0.7],...
            'ForegroundColor',[0 0 0],...
            'Position', pos.SLD);
        
        % Create select and abort push buttons
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Select',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        AbortButton = uicontrol(f,'Style', 'pushbutton', 'String', 'No Tinnitus',...
            'Position', pos.AbortBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'backgroundcolor',[1 0 0],...
            'ForegroundColor',[1 1 1],...
            'Callback', @abort);
        
        % 		Question Text
        txt = uicontrol(f,'Style','text',...
            'Position',pos.QTXT,...
            'String','How loud is your tinnitus?',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 72,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        %     0, ‘extremely weak’; 30, ‘moderate’; 50, ‘strong’; 70, ‘very strong’; and 100, ‘extremely strong’ [15].
        %         pos.ExWeakTXT
        
        ExWeakTXT = uicontrol(f,'Style','text',...
            'Position',pos.ExWeakTXT,...
            'String','extremely soft',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        moderateTXT = uicontrol(f,'Style','text',...
            'Position',pos.ModTXT,...
            'String','moderate',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        strongTXT = uicontrol(f,'Style','text',...
            'Position',pos.StrongTXT,...
            'String','strong',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        VstrongTXT = uicontrol(f,'Style','text',...
            'Position',pos.VStrongTXT,...
            'String','very strong',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        ExstrongTXT = uicontrol(f,'Style','text',...
            'Position',pos.ExStrongTXT,...
            'String','extremely loud',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        %     % Create slider value display
        val = get (sld,'Value');
        DispVal = uicontrol(f,'Style','text',...
            'Position',pos.VALTXT,...
            'String',val,...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 20,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        
        SlideLeft = uicontrol(f, 'Style', 'slider','Callback',@decrease);
        SlideRight = uicontrol(f, 'Style', 'slider','Callback',@increase);
        
        UpdateVal = uicontrol(f, 'Style', 'text','Callback',@upval);
        saveas(gcf,'PerceivedLoudness.png')
        function upval(src, e)
            val = get (sld, 'Value');
            set (DispVal, 'String', val);
            
        end
        function increase(src, e)
            
            val = get (sld,'Value');
            max = get (sld, 'Max');
            step = 1;
            
            if val < max
                valupdate = val + step;
                set (sld, 'Value', valupdate);
                %         disp (val)
            end
            
        end
        function decrease(src, e)
            
            val = get (sld,'Value');
            min = get (sld, 'Min');
            step = 1;
            
            if val > min
                valupdate = val - step;
                set (sld, 'Value', valupdate);
                %         disp (val)
            end
            
        end
        function select(src, e)
            
            val = get (sld, 'Value');
            setappdata(df,'Response',val);
            setappdata(df,'TinStat',1);
            
            close (f);
            
        end
        function abort(src, e)
            
            val = get (sld, 'Value');
            setappdata(df,'Response',val);
            setappdata(df,'TinStat',0);
            
            close (f);
            
            
        end
        function keyPress(src, e)
            
            switch e.Key
                case 'leftarrow'
                    decrease(SlideLeft, []);
                    upval(UpdateVal, []);
                    
                case 'rightarrow'
                    increase(SlideRight, []);
                    upval(UpdateVal, []);
                    
                case 'space'
                    select(SelectButton, []);
                    
                case 'a'
                    abort(AbortButton, []);
                    
            end
            
        end
        
        function pos = genGuiLayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtQSizex = scnsize(3).*0.75;
            txtQSizey = scnsize(4).*0.4;
            
            
            pos.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            pos.SLD = [scnsize(3)./2-((scnsize(3)./2)/2),...
                scnsize(4)./3,...
                scnsize(3)./2,...
                scnsize(4)./5];
            
            txtPosx = pos.SLD(1);
            txtPosy = pos.SLD(2) + (scnsize(4).*0.13);
            
            txtSizex = pos.SLD(3).*0.2;
            txtSizey = pos.SLD(4).*0.1;
            
            pos.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.25) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            pos.AbortBTN = [(scnsize(3).*0.75) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            pos.QTXT = [(scnsize(3).*0.5) - (txtQSizex./2),...
                (scnsize(4).*0.65) - (txtQSizey./2),...
                txtQSizex,...
                txtQSizey];
            
            %     0, ‘extremely weak’; 30, ‘moderate’; 50, ‘strong’; 70, ‘very strong’; and 100, ‘extremely strong’ [15].
            pos.ExWeakTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.08),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.ModTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.32),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.StrongTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.5),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.VStrongTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.69),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.ExStrongTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.91),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.VALTXT = [(scnsize(3).*0.25) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
        end
        % end of functions within GUI
    end
% end of functions within perceived loudness
end
function data = localize (RP, expData)
disp('Tinnitus characterisation ... Localize')
done = 0;
data.Localize = NaN (1);
data.InduceTime = 0;
while (done == 0)
    %% Start TDT playback HERE
    LocalizePlayer (RP, expData);
    
    %% Load GUI for participant response and parameter updates
    [f,df] = LocalizeGUI (RP, expData);
    uiwait(f)
    TinStat = getappdata(df,'TinStat');
    
    data.Localize = getappdata(df,'Response');
    
    if (TinStat == 0)
        %% run inducing stimuli
        disp ('Tinnitus no longer present... Inducing tinnitus')
        TinIn = 0;
        while (TinIn == 0)
            %             duration = 10; % Time in minutes to present tinnitus inducing stimulus
            TinntusInduce (RP, expData, expData.repDuration);
            data.InduceTime = data.InduceTime + expData.repDuration;
            
            %% Check if tinnitus has been induced
            % Load yes/no GUI
            TinIn = checkTinStat (expData);
            %             if (TinIn == 1)
            %             break
            %             end
        end
        
    elseif (TinStat == 1)
        disp ('Localization complete')
        done = 1;
    end
    
end

    function [f, df] = LocalizeGUI (RP, expData)
        
        %% Data transfer figure
        % Create a figure to save data to
        df = figure ('Visible','off',...
            'Menu', 'None');
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genLocalizeGUILayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        
        setappdata(df,'TinStat',1);
        
        
        %% Create GUI elements
        
        % Create slider
        sld = uicontrol(f,'Style', 'slider',...
            'Min',-1,'Max',1,'Value',0,...
            'backgroundcolor',[0.7 0.7 0.7],...
            'ForegroundColor',[0 0 0],...
            'Visible','off',...
            'Position', pos.SLD);
        
        % Create select and abort push buttons
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Select',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        AbortButton = uicontrol(f,'Style', 'pushbutton', 'String', 'No Tinnitus',...
            'Position', pos.AbortBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'backgroundcolor',[1 0 0],...
            'ForegroundColor',[1 1 1],...
            'Callback', @abort);
        
        % 		Question Text
        txt = uicontrol(f,'Style','text',...
            'Position',pos.QTXT,...
            'String','Which ear is your tinnitus coming from?',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 72,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        LEFTtxt = uicontrol(f,'Style','text',...
            'Position',pos.LeftTXT,...
            'String','Left',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'FontSize', 65,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        BOTHtxt = uicontrol(f,'Style','text',...
            'Position',pos.BothTXT,...
            'String','Both',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor','cyan',...
            'FontSize', 65,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        RIGHTtxt = uicontrol(f,'Style','text',...
            'Position',pos.RightTXT,...
            'String','Right',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'FontSize', 65,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        
        SlideLeft = uicontrol(f, 'Style', 'slider','Callback',@decrease);
        SlideRight = uicontrol(f, 'Style', 'slider','Callback',@increase);
        saveas(gcf,'LocalizeGUI.png')
        %% Play sound
        RP.SoftTrg(1);
        
        function increase(src, e)
            
            val = get (sld,'Value');
            max = get (sld, 'Max');
            step = 1;
            
            if val < max
                
                valupdate = val + step;
                
                CurLoc (valupdate);
                set (sld, 'Value', valupdate);
                
            end
            
        end
        function decrease(src, e)
            val = get (sld,'Value');
            min = get (sld, 'Min');
            step = 1;
            
            if val > min
                
                valupdate = val - step;
                
                CurLoc (valupdate);
                set (sld, 'Value', valupdate);
                
            end
            
        end
        
        function select(src, e)
            val = get (sld, 'Value');
            setappdata(df,'Response',val);
            setappdata(df,'TinStat',1);
            
            % Stop playing
            RP.SoftTrg(2);
            
            close (f);
        end
        function abort(src, e)
            val = get (sld, 'Value');
            setappdata(df,'Response',val);
            setappdata(df,'TinStat',0);
            
            % Stop playing
            RP.SoftTrg(2);
            
            close (f);
        end
        function keyPress(src, e)
            
            
            switch e.Key
                case 'leftarrow'
                    decrease(SlideLeft, []);
                    
                    
                case 'rightarrow'
                    increase(SlideRight, []);
                    
                    
                case 'space'
                    select(SelectButton, []);
                    
                case 'a'
                    abort(AbortButton, []);
                    
            end
            
        end
        
        function CurLoc (val)
            
            switch val
                case -1
                    disp(-1)
                    set (LEFTtxt, 'ForegroundColor', 'cyan');
                    set (BOTHtxt, 'ForegroundColor', 'k');
                    set (RIGHTtxt, 'ForegroundColor', 'k');
                    
                    RP.SoftTrg(2);
                    pause(0.2)
                    RP.SetTagVal('StereoChannel',1);
                    RP.SetTagVal('OutputChannel',2);
                    %                 pause(0.1)
                    RP.SoftTrg(1);
                    
                case 0
                    disp(0)
                    set (LEFTtxt, 'ForegroundColor', 'k');
                    set (BOTHtxt, 'ForegroundColor', 'cyan');
                    set (RIGHTtxt, 'ForegroundColor', 'k');
                    
                    RP.SoftTrg(2);
                    pause(0.2)
                    RP.SetTagVal('StereoChannel',0);
                    RP.SetTagVal('OutputChannel',1);
                    %                 pause(0.1)
                    RP.SoftTrg(1);
                    
                case 1
                    disp(1)
                    set (LEFTtxt, 'ForegroundColor', 'k');
                    set (BOTHtxt, 'ForegroundColor', 'k');
                    set (RIGHTtxt, 'ForegroundColor', 'cyan');
                    
                    RP.SoftTrg(2);
                    pause(0.2)
                    RP.SetTagVal('StereoChannel',2);
                    RP.SetTagVal('OutputChannel',2);
                    %                 pause(0.1)
                    RP.SoftTrg(1);
            end
            
        end
        
        
        function data = genLocalizeGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtQSizex = scnsize(3).*0.75;
            txtQSizey = scnsize(4).*0.4;
            txtSizex = scnsize(3).*0.2;
            txtSizey = scnsize(4).*0.2;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.SLD = [scnsize(3)./2-((scnsize(3)./2)/2),...
                scnsize(4)./3,...
                scnsize(3)./2,...
                scnsize(4)./5];
            
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.25) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            data.AbortBTN = [(scnsize(3).*0.75) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            data.QTXT = [(scnsize(3).*0.5) - (txtQSizex./2),...
                (scnsize(4).*0.65) - (txtQSizey./2),...
                txtQSizex,...
                txtQSizey];
            
            data.LeftTXT = [(scnsize(3).*0.25) - (txtSizex./2),...
                (scnsize(4).*0.4) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.BothTXT = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.4) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.RightTXT = [(scnsize(3).*0.75) - (txtSizex./2),...
                (scnsize(4).*0.4) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
        end
        % end of functions within GUI
    end
    function LocalizePlayer (RP, expData)
        
        stim_dB = 70;
        stim_amp=10^((stim_dB-expData.NomLev-expData.HB7Gain)/20);
        Freq = 500;
        RP.SetTagVal('ToneFreq',Freq)
        RP.SetTagVal('SourceChannel',0);
        RP.SetTagVal('Amp',stim_amp);
        RP.SetTagVal('StereoChannel',0);
        RP.SetTagVal('OutputChannel',1);
        
    end

end
function data = likeness (RP, expData, Spectral, LDnessRatings)
disp('Tinnitus characterisation ... Likeness')
%% n = number of frequencies to test
% nr = number of times to repeat each frequency presentation
n = 8 ;
nr = 3 ;
%% Generate test order
% output
% data.iTrial
% data.TinTestFreqOrder
test = genExpPar (n, nr);

%% while loop until all frequencies presented
done = 0;
q = 1;
index = zeros (n);
data.TinRes = NaN (n,nr);
data.InduceTime = 0;

while (done == 0)
    for i = (q:test.iTrial)
        
        c = test.TinTestFreqOrder(i:end);
        cS = c(1);
        
        switch (cS)
            case 1
                freq = 250;
                loudnessCur = LDnessRatings(1,1);
                index (cS) = index (cS) + 1;
            case 2
                freq = 500;
                loudnessCur = LDnessRatings(2,1);
                index (cS) = index (cS) + 1;
            case 3
                freq = 1000;
                loudnessCur = LDnessRatings(3,1);
                index (cS) = index (cS) + 1;
            case 4
                freq = 2000;
                loudnessCur = LDnessRatings(4,1);
                index (cS) = index (cS) + 1;
            case 5
                freq = 4000;
                loudnessCur = LDnessRatings(5,1);
                index (cS) = index (cS) + 1;
            case 6
                freq = 6000;
                loudnessCur = LDnessRatings(6,1);
                index (cS) = index (cS) + 1;
            case 7
                freq = 8000;
                loudnessCur = LDnessRatings(7,1);
                index (cS) = index (cS) + 1;
            case 8
                freq = 12000;
                loudnessCur = LDnessRatings(8,1);
                index (cS) = index (cS) + 1;
        end
        
        disp(freq)
        
        %% Start TDT playback
        LikenessPlayer (RP, expData, freq, Spectral, loudnessCur);
        
        %% Load GUI for participant response and parameter updates
        [f,df] = LikenessGUI  (expData);
        uiwait(f)
        TinStat = getappdata(df,'TinStat');
        
        % Save to first nan col for the row for that frequency(cS)
        TinResponse = getappdata(df,'Response');
        
        col =  index (cS);
        data.TinRes(cS, col) = TinResponse;
        
        if (TinStat == 0)
            %% run inducing stimuli
            disp ('Tinnitus no longer present... Inducing tinnitus')
            q = i;
            TinIn = 0;
            while (TinIn == 0)
                %             duration = 10; % Time in minutes to present tinnitus inducing stimulus
                TinntusInduce (RP, expData, expData.repDuration);
                data.InduceTime = data.InduceTime + expData.repDuration;
                %% Check if tinnitus has been induced
                % Load yes/no GUI
                TinIn = checkTinStat (expData);
                if (TinIn == 1)
                    break
                end
            end
        elseif (TinStat == 1)
            disp ('Tinnitus still present... continue characterisation')
        end
        
    end
    if (i == test.iTrial && TinStat == 1 )
        done = 1;
        disp ('Characterisation complete')
    end
    
end

close (df)

    function [f,df] = LikenessGUI  (expData)
        
        %% Data transfer figure
        % Create a figure to save data to
        df = figure ('Visible','off',...
            'Menu', 'None');
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genGuiLayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        setappdata(df,'TinStat',1);
        
        %% Create GUI elements
        % Create slider
        sld = uicontrol(f,'Style', 'slider',...
            'Min',0,'Max',100,'Value',0,...
            'backgroundcolor',[0.7 0.7 0.7],...
            'ForegroundColor',[0 0 0],...
            'Position', pos.SLD);
        
        % Create select and abort push buttons
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Select',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        AbortButton = uicontrol(f,'Style', 'pushbutton', 'String', 'No Tinnitus',...
            'Position', pos.AbortBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'backgroundcolor',[1 0 0],...
            'ForegroundColor',[1 1 1],...
            'Callback', @abort);
        
        % 		Question Text
        txt = uicontrol(f,'Style','text',...
            'Position',pos.QTXT,...
            'String','How similar is the PITCH of the sound to your tinnitus?',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 72,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        %     0, ‘extremely weak’; 30, ‘moderate’; 50, ‘strong’; 70, ‘very strong’; and 100, ‘extremely strong’ [15].
        
        ExWeakTXT = uicontrol(f,'Style','text',...
            'Position',pos.ExWeakTXT,...
            'String','Not at all',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        moderateTXT = uicontrol(f,'Style','text',...
            'Position',pos.ModTXT,...
            'String','Not very',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        strongTXT = uicontrol(f,'Style','text',...
            'Position',pos.StrongTXT,...
            'String','Somewhat similar',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        VstrongTXT = uicontrol(f,'Style','text',...
            'Position',pos.VStrongTXT,...
            'String','Very similar',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        ExstrongTXT = uicontrol(f,'Style','text',...
            'Position',pos.ExStrongTXT,...
            'String','Identical',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 12,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        %     % Create slider value display
        val = get (sld,'Value');
        DispVal = uicontrol(f,'Style','text',...
            'Position',pos.VALTXT,...
            'String',val,...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 20,...
            'FontUnits', 'normalized',...
            'FontWeight','normal');
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        
        SlideLeft = uicontrol(f, 'Style', 'slider','Callback',@decrease);
        SlideRight = uicontrol(f, 'Style', 'slider','Callback',@increase);
        
        UpdateVal = uicontrol(f, 'Style', 'text','Callback',@upval);
        saveas(gcf,'Likeness.png')
        %% Play sound
        RP.SoftTrg(1);
        
        function upval(src, e)
            val = get (sld, 'Value');
            set (DispVal, 'String', val);
            
        end
        function increase(src, e)
            
            val = get (sld,'Value');
            max = get (sld, 'Max');
            step = 1;
            
            if val < max
                valupdate = val + step;
                set (sld, 'Value', valupdate);
                %         disp (val)
            end
            
        end
        function decrease(src, e)
            
            val = get (sld,'Value');
            min = get (sld, 'Min');
            step = 1;
            
            if val > min
                valupdate = val - step;
                set (sld, 'Value', valupdate);
                %         disp (val)
            end
            
        end
        function select(src, e)
            
            val = get (sld, 'Value');
            setappdata(df,'Response',val);
            setappdata(df,'TinStat',1);
            
            % Stop playing
            RP.SoftTrg(2);
            pause(0.3)
            close (f);
            pause(1.5)
            
        end
        function abort(src, e)
            
            val = get (sld, 'Value');
            %             setappdata(df,'Response',val);
            setappdata(df,'TinStat',0);
            
            % Stop playing
            RP.SoftTrg(2);
            pause(0.3)
            close (f);
            pause(1.5)
            
        end
        function keyPress(src, e)
            
            switch e.Key
                case 'leftarrow'
                    decrease(SlideLeft, []);
                    upval(UpdateVal, []);
                    
                case 'rightarrow'
                    increase(SlideRight, []);
                    upval(UpdateVal, []);
                    
                case 'space'
                    select(SelectButton, []);
                    
                case 'a'
                    abort(AbortButton, []);
                    
            end
            
        end
        
        function pos = genGuiLayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtQSizex = scnsize(3).*0.75;
            txtQSizey = scnsize(4).*0.4;
            
            
            pos.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            pos.SLD = [scnsize(3)./2-((scnsize(3)./2)/2),...
                scnsize(4)./3,...
                scnsize(3)./2,...
                scnsize(4)./5];
            
            txtPosx = pos.SLD(1);
            txtPosy = pos.SLD(2) + (scnsize(4).*0.13);
            
            txtSizex = pos.SLD(3).*0.2;
            txtSizey = pos.SLD(4).*0.1;
            
            pos.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.25) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            pos.AbortBTN = [(scnsize(3).*0.75) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            pos.QTXT = [(scnsize(3).*0.5) - (txtQSizex./2),...
                (scnsize(4).*0.65) - (txtQSizey./2),...
                txtQSizex,...
                txtQSizey];
            
            %     0, ‘extremely weak’; 30, ‘moderate’; 50, ‘strong’; 70, ‘very strong’; and 100, ‘extremely strong’ [15].
            pos.ExWeakTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.08),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.ModTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.32),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.StrongTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.5),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.VStrongTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.69),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.ExStrongTXT = [(txtPosx - (txtSizex./2))+(pos.SLD(3).*0.91),...
                (txtPosy - (txtSizey./2)),...
                txtSizex,...
                txtSizey];
            
            pos.VALTXT = [(scnsize(3).*0.25) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
        end
        % end of functions within GUI
    end

    function LikenessPlayer (RP, expData, freq, Spectral, loudnessCur)
        
        %     stim_dB = loudnessCur;
        
        stim_amp=10^((loudnessCur-expData.NomLev-expData.HB7Gain)/20);
        
        switch Spectral
            case -1
                Switch = 0;
                RP.SetTagVal('ToneFreq',freq);
            case 0
                BW = ((freq/100).*2);
                LPfreq = freq + BW;
                HPfreq = freq - BW;
                RP.SetTagVal('LPFilterFreq',LPfreq);
                RP.SetTagVal('HPFilterFreq',HPfreq);
                FilterComp = 10*log10(BW/(expData.SF/2));
                
                %             RP.SetTagVal('FilterComp',FilterComp);
                
                Switch = 1;
            case 1
                BW = ((freq/100).*10);
                LPfreq = freq + BW;
                HPfreq = freq - BW;
                RP.SetTagVal('LPFilterFreq',LPfreq);
                RP.SetTagVal('HPFilterFreq',HPfreq);
                
                FilterComp = 10*log10(BW/(expData.SF/2));
                
                %             RP.SetTagVal('FilterComp',FilterComp);
                
                Switch = 1;
        end
        
        RP.SetTagVal('SourceChannel',Switch);
        RP.SetTagVal('Amp',stim_amp);
        RP.SetTagVal('StereoChannel',0);
        RP.SetTagVal('OutputChannel',1)
        
    end
% end of functions within likeness function
end
%% Instructions
function StartInstructions (expData)
disp('Start Instructions')

%% Load GUI for participant response and parameter updates
[f] = StartInstructionsGUI  (expData);
uiwait(f)


    function [f] = StartInstructionsGUI (expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genLocalizeGUILayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        %% Create GUI elements
        
        
        % Create select and abort push buttons
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Select',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        % 		Question Text
        txt1 = uicontrol(f,'Style','text',...
            'Position',pos.TXT1,...
            'String','To start, we will record your resting state EEG.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        txt2 = uicontrol(f,'Style','text',...
            'Position',pos.TXT2,...
            'String','A fixation cross will appear on screen. Please focus on the cross and stay as still as possible.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        %
        txt3 = uicontrol(f,'Style','text',...
            'Position',pos.TXT3,...
            'String','When ready, press the dial to continue.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        
        saveas(gcf,'Startinst.png')
        function select(src, e)
            close (f);
        end
        
        function keyPress(src, e)
            
            switch e.Key
                case 'space'
                    select(SelectButton, []);
                    
            end
            
        end
        
        function data = genLocalizeGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.75;
            txtSizey = scnsize(4).*0.3;
            txtSizey1 = scnsize(4).*0.5;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            data.TXT1 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.7) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.TXT2 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.4) - (txtSizey1./2),...
                txtSizex,...
                txtSizey1];
            
            data.TXT3 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.2) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
        end
        % end of functions within GUI
    end

end
function InduceInstructions (expData)
disp('Induce Instructions')

%% Load GUI for participant response and parameter updates
[f] = InduceInstructionsGUI  (expData);
uiwait(f)


    function [f] = InduceInstructionsGUI (expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genLocalizeGUILayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        %% Create GUI elements
        
        
        % Create select and abort push buttons
        SelectButton = uicontrol(f,'Style', 'pushbutton', 'String', 'Select',...
            'Position', pos.SelectBTN,...
            'FontSize', 25,...
            'FontUnits', 'normalized',...
            'FontWeight','normal',...
            'Callback', @select);
        
        % 		Question Text
        txt1 = uicontrol(f,'Style','text',...
            'Position',pos.TXT1,...
            'String','We will now present the auditory stimuli.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        txt2 = uicontrol(f,'Style','text',...
            'Position',pos.TXT2,...
            'String','Once the presentation is complate, a fixation cross will appear on screen. Please focus on the cross and stay as still as possible.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        %
        txt3 = uicontrol(f,'Style','text',...
            'Position',pos.TXT3,...
            'String','When ready, press the dial to continue.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        
        %% Create key control of GUI elements
        
        set (f,'KeyPressFcn', @keyPress)
        
        saveas(gcf,'InduceInst.png')
        function select(src, e)
            close (f);
        end
        
        function keyPress(src, e)
            
            switch e.Key
                case 'space'
                    select(SelectButton, []);
                    
            end
            
        end
        
        function data = genLocalizeGUILayout (f, expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            butSizeX = scnsize(3)./10;
            butSizey = scnsize(4)./20;
            txtSizex = scnsize(3).*0.75;
            txtSizey = scnsize(4).*0.3;
            txtSizey1 = scnsize(4).*0.5;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            
            data.SelectBTN = [(scnsize(3).*0.5) - (butSizeX./2),...
                (scnsize(4).*0.15) - (butSizey./2),...
                butSizeX,...
                butSizey];
            
            data.TXT1 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.7) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.TXT2 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.4) - (txtSizey1./2),...
                txtSizex,...
                txtSizey1];
            
            data.TXT3 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.1) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
        end
        % end of functions within GUI
    end

end

function ThankYouScreen (expData)
disp('Experiment Complete')

%% Load GUI for participant response and parameter updates
[f] = ThankYouScreenGUI  (expData);
prompt = 'Hit Enter to Close figure';
str = input(prompt);
close (f)

    function [f] = ThankYouScreenGUI (expData)
        
        % Create a figure for present
        f = figure ('Visible','off',...
            'Menu', 'None', ...
            'Color', [0.5, 0.5, 0.5]);
        
        pos = genLocalizeGUILayout (f, expData);
        
        set(f,'OuterPosition',pos.FIG)
        
        %% Create GUI elements
        
        % 	Text
        txt1 = uicontrol(f,'Style','text',...
            'Position',pos.TXT1,...
            'String','Thank you for taking part in this study.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        txt2 = uicontrol(f,'Style','text',...
            'Position',pos.TXT2,...
            'String','Please wait for the experimenter to enter the room.',...
            'backgroundcolor',[0.5, 0.5, 0.5],...
            'ForegroundColor',[1 1 1],...
            'FontSize', 60,...
            'FontUnits', 'normalized',...
            'FontWeight','bold');
        
        
        % Make figure visble after adding all components
        set (f, 'Visible', 'on');
        saveas(gcf,'Thanks.png')
        
        function data = genLocalizeGUILayout (expData)
            
            scnsize = expData.scnsize;
            edge = expData.edge;
            
            txtSizex = scnsize(3).*0.95;
            txtSizey = scnsize(4).*0.3;
            txtSizey1 = scnsize(4).*0.5;
            
            data.FIG = [scnsize(1),...
                scnsize(2),...
                scnsize(3) - edge,...
                scnsize(4)];
            
            data.TXT1 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.7) - (txtSizey./2),...
                txtSizex,...
                txtSizey];
            
            data.TXT2 = [(scnsize(3).*0.5) - (txtSizex./2),...
                (scnsize(4).*0.45) - (txtSizey1./2),...
                txtSizex,...
                txtSizey1];
            
            
        end
        % end of functions within GUI
    end

end
