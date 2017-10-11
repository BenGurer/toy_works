function [stimInfo, glmInfo] = sHL_setupStudyParams
    %
    %   usage: sHL_setupStudyParams
    %      by: Ben Gurer
    %    date: 11/10/2017
    % purpose: setup study info and params for hearing loss simulation study
    %   input: n/a
    %  output: stimInfo, stimulus information; glmInfo, glm analysis information
    %
    
%% Get stimulus properties
% get stimulus frequencies
lowFreqkHz = 0.1;
highFreqkHz = 8;
nStim = 32;

[stimInfo.stimFreqs, stimInfo.stimFreqs_bin, stimInfo.stimFreqs_mv] = convertStimIDtoFrequency(lowFreqkHz,highFreqkHz,nStim);

%% get stimulus senssation level
[stimLevel_SL, maskingLevel] = calStimulusSensationLevel(stimInfo.stimFreqs);

% bin sensation level
binSize = 4;
stimLevel_SL_bin = zeros(1,length(stimLevel_SL)/binSize);
c = 1;
for i = 1:length(stimLevel_SL)/binSize
    stimLevel_SL_bin(i) = mean(stimLevel_SL(c:c + (binSize-1)));
    c = c + binSize;
end

% moving average sensation level
nBins = 8;
windowAvSize = length(stimLevel_SL)/nBins;

if isreal(windowAvSize) && rem(windowAvSize,1)==0
    loopLength = length(stimLevel_SL) - windowAvSize;
    stimLevel_SL_mv = zeros(1,loopLength);
    for i = 1:loopLength
        stimLevel_SL_mv(i) = mean(stimLevel_SL(i:i+windowAvSize-1));
    end
else
    error('Moving average window not an integer')
end

% save in structure
stimInfo.stimLevel_SL = stimLevel_SL;
stimInfo.stimLevel_SL_bin = stimLevel_SL_bin;
stimInfo.stimLevel_SL_mv = stimLevel_SL_mv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup GLM analysis
% save in glmInfo structure
glmInfo.hrfModel = {'hrfBoxcar', 'hrfDoubleGamma'};
glmInfo.groupNames = {'ConcatenationHLsim', 'ConcatenationNH'};
glmInfo.nScans = 4;
glmInfo.nStim = [32, 8];
glmInfo.analysisNames_Scans = cell(1,nScans);
for iStim = 1:glmInfo.nScans
for iScan = 1:nScans
    for iHRF = 1:length(glmInfo.hrfModel)
        glmInfo.analysisNames_Scans{iScan} = ['glm_' glmInfo.hrfModel{iHRF} '_nCons_' mat2str(glmInfo.nStim(iStim)) '_Scan_' mat2str(iScan)];
    end
end
end

glmInfo.analysisNames_Groups = cell(1,length(glmInfo.groupNames));
for iGroup = 1:length(glmInfo.groupNames)
    for iHRF = 1:length(glmInfo.hrfModel)
        glmInfo.analysisNames_Groups = ['glm_' glmInfo.hrfModel{iHRF}];
    end
end