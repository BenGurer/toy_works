function overlayWeighted = weightBetaOverlay(overlay,fit)

if fit == 0
    weightingType = 'SL';
else
    weightingType = 'BOLD';
end

weight = stimWeighting(weightingType,fit);

overlayWeighted =  overlay./ weight;

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