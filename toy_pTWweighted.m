

% use a drop down menu to choose weighting type
% use tick box to fit mod
% if nStim == stimInfo.sizes(1)
    x = stimInfo.stimNames.bin;
% elseif nStim == stimInfo.sizes(2)
%     x = stimInfo.stimNames.mv;
% elseif nStim == stimInfo.sizes(3)
% x = stimInfo.stimNames.all;
% end

stimulusLevel_dbSPL = 75;
maskingLevel_dbSPL = 25;
threshold_sHL_dBSLP = funSimulateHearingLoss(x);
masking_Baseline = maskingLevel_dbSPL*ones(size(x));
masking_dbSPL =  max(threshold_sHL_dBSLP,masking_Baseline);
stimulusLevel_dbSL = stimulusLevel_dbSPL-masking_dbSPL;

        stimulusWeighting = (stimulusLevel_dbSL)/max(stimulusLevel_dbSL);
  data2test = data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW;
  
%%
    figure
for iCon = 1:length(data2test) 
for i = 1:length(data2test{iCon})

% plot(data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_pTW{iCon}(i,:))
subplot(2,4,i)
plot(data2test{iCon}(i,:))
    hold on
    legend('show')
end
end

for i = 1:length(data2test{2})

% plot(data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_pTW{iCon}(i,:))
subplot(2,4,i)
plot(data2test{2}(i,:)./stimulusWeighting)
    hold on
end


%%

    figure
for iCon = 1:length(data2test) 
for i = 1:length(data2test{iCon})

% plot(data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_pTW{iCon}(i,:))
subplot(2,4,i)
plot(data2test{iCon}(i,:)./max(max(data2test{iCon})))
    hold on
    legend('show')
end
end
for i = 1:length(data2test{2})

% plot(data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_pTW{iCon}(i,:))
subplot(2,4,i)
plot(data2test{2}(i,:)./max(max(data2test{iCon}))./stimulusWeighting)
    hold on
end

%%
figure
for iCon = 1:length(data2test)
    subplot(1,3,iCon)
surf(data2test{iCon})
end
subplot(1,3,3)
surf(data2test{2}./stimulusWeighting)