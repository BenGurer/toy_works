function thisView = script_glmOverlayAnalysis(thisView,overlayNumbers,namePrefix,fit,conditionNames)

%% Correction of beta weights

weight = weightBetaOverlay(fit);
overlayNames = cell(1,length(overlayNumbers));
for iOverlay = 1:length(overlayNumbers)
    [thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlayNumbers(iOverlay))]);
    params.combineFunction='User Defined';
    params.customCombineFunction = 'rdivide'; % (overlay,weightingType,fit)
    params.combinationMode = 'Apply function to each overlay';
    params.outputName= namePrefix;
    eval(['params.additionalArgs = ' num2str(weight(iOverlay)) ';']);
    % params.baseSpace = 1;
    [thisView,params] = combineTransformOverlays(thisView,params);
    overlayNames{iOverlay} = [namePrefix '(' conditionNames{iOverlay} ',' num2str(weight(iOverlay)) ')'];
end

% get overlays to perform analysis on
overlayNumbers_weighted = viewGet(thisView,'overlayNum',overlayNames);

%Tonotopy analysis
% Index max
[thisView,params] = combineTransformOverlays(thisView,[],'justGetParams=1','defaultParams=1',['overlayList=' mat2str(overlayNumbers_weighted)]);
params.combineFunction='indexMax';
params.nOutputOverlays=2;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-1);

% Weighted mean and corrected weighted mean
params.combineFunction='weightedMeanStd';
params.nOutputOverlays=4;
[thisView,params] = combineTransformOverlays(thisView,params);
curOverlay=viewGet(thisView,'curOverlay');
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-3);
thisView = viewSet(thisView,'overlaycolorrange',[0 32],curOverlay-2);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay-1);
thisView = viewSet(thisView,'overlaycolorrange',[0 40],curOverlay);

end

function weight = weightBetaOverlay(fit)

if fit == 0
    weightingType = 'SL';
else
    weightingType = 'BOLD';
end

weight = stimWeighting(weightingType,fit);

end


function stimulusWeighting = stimWeighting(weightingType,fit)


[stimInfo, ~, ~, ~, ~] = sHL_setupStudyParams;
%
% use a drop down menu to choose weighting type
% use tick box to fit mod
% if nStim == stimInfo.sizes(1)
%     x = stimInfo.stimNames.bin;
% elseif nStim == stimInfo.sizes(2)
%     x = stimInfo.stimNames.mv;
% elseif nStim == stimInfo.sizes(3)
x = stimInfo.stimNames.all;
% end

stimulusLevel_dbSPL = 75;
maskingLevel_dbSPL = 25;
threshold_sHL_dBSLP = funSimulateHearingLoss(x);
masking_Baseline = maskingLevel_dbSPL*ones(size(x));
masking_dbSPL =  max(threshold_sHL_dBSLP,masking_Baseline);
stimulusLevel_dbSL = stimulusLevel_dbSPL-masking_dbSPL;

switch weightingType
    case 'SL'
        
        % %     Convert to pressure
        % %   Lp(dB SPL) = 20 log10 p/p0
        % %   p0 = 0.00002 pa
        % %   p(Pa) = p0 .10.^ Lp(dB SPL)/20
        % stimWeightingPressure = 0.00002 .* 10.^(stimulusLevel_dbSL/20);
        % %       stimWeightingIntensity = 10.^-12 .* 10.^((20 .* log10(stimWeightingPressure/10.^-5))/10);
        % stimWeightingIntensity = (stimWeightingPressure.^2) / 400;
        % stimulusLevel_dbSL = stimWeightingIntensity;
        % %     stimWeighting = 10.^(stimWeighting/10);
        % stimulusLevel_dbSL = stimulusLevel_dbSL/max(stimulusLevel_dbSL);
        % %     stimWeighting = (stimLevel-threshEvel)/stimSLlevel;
        
        % % normalise by max value
        stimulusWeighting = (stimulusLevel_dbSL)/max(stimulusLevel_dbSL);
        
    case 'BOLD'
        %         R = m.SL + b;
        m = fit(1);
        b = fit(2);
        %         m = 0.0174;
        %         b = -0.1176;
        R = @(SL) m.*SL + b;
        stimulusWeighting = R(stimulusLevel_dbSL);
        
    case 'fit'
        
        stimulusWeighting = (stimulusLevel_dbSL)/max(stimulusLevel_dbSL);
        stimulusWeighting = stimulusWeighting.^compression;
end
end

