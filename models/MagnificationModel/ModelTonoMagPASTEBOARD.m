function ModelTonoMagPASTEBOARD
if fESTpCFPrunedSel == 1;
    %% pCF lim StimFreqs
    % create figure plotting estimated pCF for each magnification domain for
    % each stimulus set with axis limits set by stimulus set min and maximum
    for iTonoMap = 1:length(TonotopicMagnification)
        fTonotopicZoomPrune(iTonoMap) = figure('Color', [1 1 1]);
        set(fTonotopicZoomPrune(iTonoMap),'OuterPosition',FigSize)
        plot (VoxelDistance,TonotopicMap(iTonoMap).pCF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2);
        axis tight, title ([TonotopicNames{iTonoMap},' - Range: ', num2str(StimLowFreq),' - ',num2str(StimHighFreq),' kHz'], 'FontSize', 20);
        xmin = find(TonotopicMap(iTonoMap).pCF > StimLowFreq,1);
        xmax = find(TonotopicMap(iTonoMap).pCF >= StimHighFreq,1);  % highest frequency to test
        pCFprune = TonotopicMap(iTonoMap).pCF(xmin+1:xmax);
        xminDis = (xmin/ nVoxels) * max(VoxelDistance);
        xmaxDis = (xmax/ nVoxels) * max(VoxelDistance);
        TonotopicMap(iTonoMap).VoxelsWithinStimLim = xmax - xmin; % need to convert from voxels to mm+
        TonotopicMap(iTonoMap).CorticalDistanceWithinStimLim = xmaxDis - xminDis;
        ylim ([StimLowFreq StimHighFreq]), xlim ([xminDis xmaxDis]);
        hold on
        SE = zeros(iStimSet,length(pCFprune));
        for iStimSet = 1:length(StimulusSetNames)
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            EstpCF = [zeros(1,xmin), EstpCF(xmin+1:xmax), zeros(1,length(EstpCF) - xmax)];
            scatter (VoxelDistance,EstpCF, 50, SpacingColour{iStimSet}, 'fill');
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            EstpCF = EstpCF(xmin+1:xmax);
            
            %% Calculate Root Mean Squared Error for voxels within stimulus frequency limits
            SE(iStimSet,:)= (log10(EstpCF) - log10(pCFprune)).^2;
            MSE(iStimSet) = mean(SE(iStimSet,:));
            TonotopicMap(iTonoMap).RMSEPruned(iStimSet) = sqrt(MSE(iStimSet));
        end
        set(gca,'YScale','log',...
            'YMinorTick','on',...
            'YTickLabel',{'0.5','1','2','4','8'},...
            'YTick',[0.5 1 2 4 8]);
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        fLeg = legend('Actual',...
            [StimulusSetNames{1},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSEPruned(1))],...
            [StimulusSetNames{2},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSEPruned(2))],...
            [StimulusSetNames{3},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSEPruned(3))],...
            [StimulusSetNames{4},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSEPruned(4))],...
            'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
    
    %% Save RMSE for voxels with in stumulus frequncy limits to excel file
    T = array2table(vertcat(TonotopicMap.RMSEPruned),'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'},...
        'RowNames',{'FreqDis' 'ERB' 'Log'});
    writetable(T,'RMSEStimPruned.csv');
end


if fRMSErange == 1;
    %% pCF lim StimFreqs
    % create figure plotting estimated pCF for each magnification domain for
    % each stimulus set with axis limits 1 to 4 kHz
    for iTonoMap = 1:length(TonotopicMagnification)
        fTonotopicZoomPrune(iTonoMap) = figure('Color', [1 1 1]);
        set(fTonotopicZoomPrune(iTonoMap),'OuterPosition',FigSize)
        plot (VoxelDistance,TonotopicMap(iTonoMap).pCF, 'color',SpacingColour{iTonoMap}, 'LineWidth', 2);
        axis tight, title ([TonotopicNames{iTonoMap},' - Range: ', num2str(rangelow),' - ',num2str(rangehigh),' kHz'], 'FontSize', 20);
        xmin = find(TonotopicMap(iTonoMap).pCF > rangelow,1);
        xmax = find(TonotopicMap(iTonoMap).pCF >= rangehigh,1);  % highest frequency to test
        pCFprune = TonotopicMap(iTonoMap).pCF(xmin+1:xmax);
        xminDis = (xmin/ nVoxels) * max(VoxelDistance);
        xmaxDis = (xmax/ nVoxels) * max(VoxelDistance);
        TonotopicMap(iTonoMap).VoxelsWithinStimLim = xmax - xmin; % need to convert from voxels to mm+
        TonotopicMap(iTonoMap).CorticalDistanceWithinStimLim = xmaxDis - xminDis;
        ylim ([StimLowFreq StimHighFreq]), xlim ([xminDis xmaxDis]);
        hold on
        SE = zeros(iStimSet,length(pCFprune));
        for iStimSet = 1:length(StimulusSetNames)
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            EstpCF = [zeros(1,xmin), EstpCF(xmin+1:xmax), zeros(1,length(EstpCF) - xmax)];
            scatter (VoxelDistance,EstpCF, 50, SpacingColour{iStimSet}, 'fill');
            EstpCF = TonotopicMap(iTonoMap).ESTpCF(iStimSet,:);
            EstpCF = EstpCF(xmin+1:xmax);
            
            %% Calculate Root Mean Squared Error for voxels within 'Range' frequency limits
            SE(iStimSet,:)= (log10(EstpCF) - log10(pCFprune)).^2;
            MSE(iStimSet) = mean(SE(iStimSet,:));
            TonotopicMap(iTonoMap).RMSERange(iStimSet) = sqrt(MSE(iStimSet));
        end
        set(gca,'YScale','log',...
            'YMinorTick','on',...
            'YTickLabel',{'0.5','1','2','4','8'},...
            'YTick',[0.5 1 2 4 8]);
        xlabel('Cortical Distance (mm)', 'FontSize', fFontSize), ylabel('pCF (kHz)', 'FontSize', fFontSize);
        fLeg = legend('Actual',...
            [StimulusSetNames{1},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSERange(1))],...
            [StimulusSetNames{2},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSERange(2))],...
            [StimulusSetNames{3},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSERange(3))],...
            [StimulusSetNames{4},' - RMSE = ',num2str(TonotopicMap(iTonoMap).RMSERange(4))],...
            'Location','SouthEast');
        set(fLeg, 'FontSize', fFontSize./2);
    end
    
    %% Save RMSE for voxels with in stumulus frequncy limits to excel file
    T = array2table(vertcat(TonotopicMap.RMSERange),'VariableNames',{'Linear' 'FreqDis' 'ERB' 'Log'},...
        'RowNames',{'FreqDis' 'ERB' 'Log'});
    writetable(T,'RMSERange.csv');
end
end