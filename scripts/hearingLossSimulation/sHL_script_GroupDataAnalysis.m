function
groupData = sHL_script_GroupDataAnalysis
% A script to load individual subject data and perform group level data analysis
%
%   usage: sHL_script_GroupDataAnalysis
%      by: Ben Gurer
%    date: 17/11/2017
% purpose: perform group level data analysis data from  hearing loss simulation study
%   input: n/a

% TO DO
% get roi av, take mean and ste
% get bold vs SL, take mean and ste - compare to mean bold vs SL
% get pCF - cal correlation, cal histogram, normalise? and plot
% plot experimental info - sHL, SL, frequency
% make flatmap figures - 2 or 3 subjects


%% get study info
[stimInfo, glmInfo, pRFInfo, Info, plotInfo] = sHL_setupStudyParams;
% stimulus info
% condition names
% nummber of subjects
nSubjects = 7;
%% get subject info
for iSub = 1:nSubjects
    % subject(iSub).Info = get_SubjectInfo_sHL(iSub);
    subjectInfo(iSub)= get_SubjectInfo_sHL(iSub);
    % Subject ID
    % flatmap names
    
    %% load data
    % move to data location then load and rename or save to subject name struct
    cd(fullfile(Info.dataDir,Info.studyDir,subjectInfo(iSub).subjectID));
    filename = [subjectInfo(iSub).subjectID '_data.mat'];
    % subject(iSub).Data = load(filename);
    subjectData(iSub)= load(filename);
end

%% analysis to compare:
% GLM
% ROI average beta weight
% ROI average tuning curves

% subjectData(1).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW

%% STIMULUS PLOTS
lineThickness = 2;
markerColors = [0 0 0];
markerSize = 8;
figure('Name', 'Sensation Level', 'color', [1 1 1]);
semilogx(stimInfo.stimNames.all,stimInfo.stimLevel_SL,...
    'LineStyle','--',...
    'LineWidth',lineThickness,...
    'Color',markerColors,...
    'Marker','o',...
    'MarkerSize',markerSize,...
    'MarkerEdgeColor',markerColors,...
    'MarkerFaceColor',markerColors)
hold on
semilogx(stimInfo.stimNames.all,repmat(50,size(stimInfo.stimNames.all)),...
    'LineStyle','--',...
    'LineWidth',lineThickness,...
    'Color',markerColors,...
    'Marker','o',...
    'MarkerSize',markerSize,...
    'MarkerEdgeColor',markerColors,...
    'MarkerFaceColor',markerColors)
xlabel('Frequency (kHz)')
ylabel('Sensation Level (db SL)')


for iSub = 1:nSubjects
    % ROIAverageTuningCurves_ConA(:,:,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{1};
    % ROIAverageTuningCurves_ConB(:,:,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{2};
    
    % ROIAverageTuningCurves(:,:,1,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{1};
    % ROIAverageTuningCurves(:,:,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{2};
    % ROIAverageTuningCurves_Norm(:,:,:,iSub) = ROIAverageTuningCurves(:,:,:,iSub) ./ max(max(max(ROIAverageTuningCurves(:,:,:,iSub))));
    
    
    ROIAverageTuningCurves(:,:,1,1,iSub) = subjectData(iSub).data.Left.LeftAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{1};
    ROIAverageTuningCurves(:,:,2,1,iSub) = subjectData(iSub).data.Left.LeftAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{2};
    ROIAverageTuningCurves_Norm(:,:,:,1,iSub) = ROIAverageTuningCurves(:,:,:,1,iSub) ./ max(max(max(ROIAverageTuningCurves(:,:,:,1,iSub))));
    
    ROIAverageTuningCurves(:,:,1,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{1};
    ROIAverageTuningCurves(:,:,2,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{2};
    ROIAverageTuningCurves_Norm(:,:,:,2,iSub) = ROIAverageTuningCurves(:,:,:,2,iSub) ./ max(max(max(ROIAverageTuningCurves(:,:,:,2,iSub))));
    
    
    % data_conA = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{1};
    % data_conB = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_8.roiAnalysis.roi_pTW{2};
    
    % data.x = [1:8,1:8];
    % data.y = [data_conA,data_conB];
    % [con{1:8}] = deal('Condition A')
    % data.con = repcell(size(data_conA))
end

% GLM pTW

%% ROI average tuning curves
% add max beta point?
% ROIAverageTuningCurves_Norm(:,:,:,5) = nan;
figx = 5;
bottom = 5;
width = 29.7;
% width = 21.0;
height =  21.0;
% height = 29.7;
xlabelnames = stimInfo.stimNames.bin;
for iSide = 1:2
    figure('Units','centimeters',...
        'Position',[figx bottom width height],...
        'color',[1 1 1],...
        'Name','pTW')
    
    minBeta = min(min(min(min(ROIAverageTuningCurves_Norm(:,:,:,iSide,:)))));
    maxBeta = max(max(max(max(ROIAverageTuningCurves_Norm(:,:,:,iSide,:)))));
    comp = 1 + (width - height)/height
    for i = 1:8
        subaxis(2,4,i,'Spacing',0,'MarginTop',0.2,'MarginBottom',0.2,'MarginLeft',0.2/comp,'MarginRight',0.2/comp)
        meanBetasA = nanmean(ROIAverageTuningCurves_Norm(i,:,1,iSide,:),5);
        meanBetasB = nanmean(ROIAverageTuningCurves_Norm(i,:,2,iSide,:),5);
        
        stdBetasA = nanstd(ROIAverageTuningCurves_Norm(i,:,1,iSide,:),0,5) ./sqrt(nSubjects);
        stdBetasB = nanstd(ROIAverageTuningCurves_Norm(i,:,2,iSide,:),0,5) ./sqrt(nSubjects);
        errorbar(meanBetasA,stdBetasA)
        hold on
        errorbar(meanBetasB,stdBetasB)

        xlim([0.5 8.5])
        ylim([-0.25 maxBeta+0.05])
        ax = gca;
        
        if i == 1 || i == 5
            ax.YTickLabelMode = 'auto';
        else
            ax.YTickLabel = [];
        end
        if i < 5
            ax.XTickLabel = [];
        else
            ax.XTickLabel = round(xlabelnames(ax.XTick),2);
        end
    end
    ax1 = axes('Position',[0.5 0.15 1 1],'Visible','off');
    axes(ax1)
    text(0,0,'Frequency (kHz)','HorizontalAlignment','center')
    
    ax2 = axes('Position',[0.1 0.5 1 1],'Visible','off');
    axes(ax2)
    text(0,0,'Beta Weight (normlised)','HorizontalAlignment','center',...
        'Rotation',90)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HERE
a = squeeze(cat(5,ROIAverageTuningCurves_Norm(:,:,:,1,:),ROIAverageTuningCurves_Norm(:,:,:,2,:)));
b = nanmean(a,4);
e = nanstd(a,0,4) ./sqrt(14);

figx = 5;
bottom = 5;
width = 29.7;
% width = 21.0;
height =  21.0;
% height = 29.7;
xlabelnames = stimInfo.stimNames.bin;
comp = 1 + (width - height)/height;
figure('Units','centimeters',...
    'Position',[figx bottom width height],...
    'color',[1 1 1],...
    'Name','pTW')
maxBeta = max(max(b));
for i = 1:8
    subaxis(2,4,i,'Spacing',0,'MarginTop',0.2,'MarginBottom',0.2,'MarginLeft',0.2/comp,'MarginRight',0.2/comp)
    p1 = plot(b(i,:,1),'color',conA_color,'LineWidth', 2)
    hold on
        errorbar(b(i,:,1),e(i,:,1),'.','color',conA_color)
    hold on
    p2 = plot(b(i,:,2),'color',conB_color,'LineWidth', 2)
    errorbar(b(i,:,2),e(i,:,2),'.','color',conB_color)
       p3= plot(i,b(i,i,1),'marker','o','color','k',...
    'MarkerFaceColor', 'k')
    xlim([0.5 8.5])
    ylim([-0.25 1])
    ax = gca;    
    ax.FontSize = 15;
    if i == 1 || i == 5
        ax.YTickLabelMode = 'auto';
    else
        ax.YTickLabel = [];
    end
    if i == 1
        legend([p1 p2 p3],'Normal Hearing', 'simulated HL', 'Max Beta','Location','northwest')
        
%         legend([p1 p2 p3],'Normal Hearing','Hearing Loss Simulation','Sensation Level','Location','southwest')
        legend('boxoff','FontSize',12)
    end
    if i < 5
        ax.XTickLabel = [];
    else
        ax.XTickLabel = round(xlabelnames(ax.XTick),1);
    end
end
ax1 = axes('Position',[0.5 0.12 1 1],'Visible','off');
axes(ax1)
text(0,0,'Frequency (kHz)','HorizontalAlignment','center','FontSize',FontSize)

ax2 = axes('Position',[0.08 0.5 1 1],'Visible','off');
axes(ax2)
text(0,0,'Beta Weight (normlised)','HorizontalAlignment','center',...
    'Rotation',90,'FontSize',FontSize)

%%% TO HERE
%%%%%%%%%%%%%
%% dB SL vs beta weight
for iSub = 1:nSubjects
    % left
    ROIAverageBetaRatio(:,1,iSub) = subjectData(iSub).data.Left.LeftAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.ratio2Plot;
    ROIAverageBetaLevel(:,1,iSub) = subjectData(iSub).data.Left.LeftAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.level2Plot;
    % ROIAverageTuningCurves_Norm(:,:,:,1,iSub) = ROIAverageTuningCurves(:,:,:,1,iSub) ./ max(max(max(ROIAverageTuningCurves(:,:,:,1,iSub))));
    % right
    ROIAverageBetaRatio(:,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.ratio2Plot;
    ROIAverageBetaLevel(:,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.level2Plot;
    % ROIAverageTuningCurves_Norm(:,:,:,2,iSub) = ROIAverageTuningCurves(:,:,:,2,iSub) ./ max(max(max(ROIAverageTuningCurves(:,:,:,2,iSub))));
    
end

figure('Units','centimeters',...
        'Position',[figx bottom width height],...
        'color',[1 1 1],...
        'Name','Bold ratio vs SL')
    
    for iSide = 1:2
        ROIAverageBetaLevel_mean = mean(ROIAverageBetaLevel(:,iSide),3);
        ROIAverageBetaRatio_mean = mean(ROIAverageBetaRatio(:,iSide),3);
        
        ROIAverageBetaLevel_ste = std(ROIAverageBetaLevel(:,iSide,:),0,3) ./sqrt(nSubjects);
        ROIAverageBetaRatio_ste = std(ROIAverageBetaRatio(:,iSide,:),0,3) ./sqrt(nSubjects);
        %% dB SL vs beta weight
        
        subaxis(1,2,iSide,'Spacing',0,'MarginTop',0.2*comp,'MarginBottom',0.2*comp,'MarginLeft',0.2/comp,'MarginRight',0.2/comp);
        errorbar( ROIAverageBetaLevel_mean , ROIAverageBetaRatio_mean , ROIAverageBetaRatio_ste,'o')
        % scatter( ROIAverageBetaLevel_mean , ROIAverageBetaRatio_mean )
        hold on
        fit = polyfit(ROIAverageBetaLevel_mean,ROIAverageBetaRatio_mean,1);
        plot(ROIAverageBetaLevel_mean,polyval(fit,ROIAverageBetaLevel_mean));
        % correlation = corrcoef([level2Plot_mv' ratio2Plot_mv]);
        % text(25,0.8,sprintf('Correlation = %.2f \n y = %.2fx + %.2f',correlation(2), fit(1), fit(2)))
        % plot(f,f,'k--')
        % errorbar(error2plot(1),error2plot(2),error2plot(3))
        xlim([14.1 50.9])
        ylim([0 1])
        % xlabel('Stimulus Sensation Level (dB SL)');
        plot(stimInfo.stimLevel_SL_mv,stimInfo.stimLevel_SL_mv./max(stimInfo.stimLevel_SL_mv),'k--')
        ax = gca;
        if iSide == 1
            ax.YTickLabelMode = 'auto';
            ylabel('Average Beta Weight Ratio (ConB / ConA)');
        else
            ax.YTickLabel = [];
        end
    end
    ax1 = axes('Position',[0.5 0.225 1 1],'Visible','off');
    axes(ax1)
    text(0,0,'Stimulus Sensation Level (dB SL)','HorizontalAlignment','center')    

%% ROI AVERAGE BETA WEIGHT
for iSub = 1:nSubjects
    
    ROIAverageBetas(:,1,1,iSub) = subjectData(iSub).data.Left.LeftAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_av{1};
    ROIAverageBetas(:,2,1,iSub) = subjectData(iSub).data.Left.LeftAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_av{2};
    ROIAverageBetas_Norm(:,:,1,iSub) = ROIAverageBetas(:,:,1,iSub) ./ max(max(max(ROIAverageBetas(:,:,1,iSub))));
    
    ROIAverageBetas(:,1,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_av{1};
    ROIAverageBetas(:,2,2,iSub) = subjectData(iSub).data.Right.RightAC.splitData.glm_hrfBoxcar_nCons_32.roiAnalysis.roi_av{2};
    ROIAverageBetas_Norm(:,:,2,iSub) = ROIAverageBetas(:,:,2,iSub) ./ max(max(max(ROIAverageBetas(:,:,2,iSub))));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ROI AVERAGE BETA WEIGHT - left and right
figure('Units','centimeters',...
    'Position',[figx bottom width height],...
    'color',[1 1 1],...
    'Name','ROI Average Beta Weight')
xlabelnamesMV = stimInfo.stimNames.mv;
for iSide = 1:2
    subaxis(1,2,iSide,'Spacing',0,'MarginTop',0.2*comp,'MarginBottom',0.2*comp,'MarginLeft',0.2/comp,'MarginRight',0.2/comp);
    for iCon = 1:2
        ROIAverageBetas_Norm_mean(:,iCon,iSide) = mean(ROIAverageBetas_Norm(:,iCon,iSide,:),4);
        ROIAverageBetas_Norm_ste(:,iCon,iSide) = std(ROIAverageBetas_Norm(:,iCon,iSide,:),0,4) ./sqrt(nSubjects);
   if iCon ==1     
        p1 = plot(1:length(ROIAverageBetas_Norm_mean(:,iCon,iSide)),ROIAverageBetas_Norm_mean(:,iCon,iSide),'LineWidth',2,'color',conA_color);
        hold on
        errorbar(1:length(ROIAverageBetas_Norm_mean(:,iCon,iSide)),ROIAverageBetas_Norm_mean(:,iCon,iSide),ROIAverageBetas_Norm_ste(:,iCon,iSide),'.','color',conA_color);
   else
        p2 = plot(1:length(ROIAverageBetas_Norm_mean(:,iCon,iSide)),ROIAverageBetas_Norm_mean(:,iCon,iSide),'LineWidth',2,'color',conB_color);
        hold on
        errorbar(1:length(ROIAverageBetas_Norm_mean(:,iCon,iSide)),ROIAverageBetas_Norm_mean(:,iCon,iSide),ROIAverageBetas_Norm_ste(:,iCon,iSide),'.','color',conB_color);
   end
        
    end% plot(mean(ROIAverageBetas_Norm(:,:,iSide),4))
    xlim([0.01, length(xlabelnamesMV)+0.99])
    ylim([0 1])
    ax = gca;
    p3 = plot(1:length(stimInfo.stimLevel_SL_mv),stimInfo.stimLevel_SL_mv/max(stimInfo.stimLevel_SL_mv),...
        'color',[ 0 0 0 ],...
        'linestyle','--',...
        'LineWidth',2);
    ax.XTickLabel = round(xlabelnamesMV(ax.XTick),2);
         title(Info.Sides{iSide},'FontSize',20)
             ax = gca;
ax.FontSize = 15;
    if iSide == 1
        ax.YTickLabelMode = 'auto';
        ylabel('Beta Weight (Normalised)','FontSize',FontSize)
        legend([p1 p2 p3],'Normal Hearing','Hearing Loss Simulation','Sensation Level','Location','southwest')
        legend('boxoff','FontSize',FontSize)
    else
        ax.YTickLabel = [];
        yyaxis right
        ylabel('Sensation Level (Normalised)','FontSize',FontSize)
        ax = gca;
        ax.YColor = [0 0 0];
        ylim([0 1])
        %         xlabel('Frequency (kHz)')
        ax1 = axes('Position',[0.5 0.22 1 1],'Visible','off');
        axes(ax1)
        text(0,0,'Frequency (kHz)','HorizontalAlignment','center','FontSize',FontSize)
    end
    %     ax.XTickLabel = round(xlabelnames(ax.XTick),2);
    %     xlabel('Frequency (kHz)')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','centimeters',...
    'Position',[figx bottom width height],...
    'color',[1 1 1],...
    'Name','Bold vs SL')
xlabelnamesMV = stimInfo.stimNames.mv;
for iSide = 1:2
    subaxis(1,2,iSide,'Spacing',0,'MarginTop',0.2*comp,'MarginBottom',0.2*comp,'MarginLeft',0.2/comp,'MarginRight',0.2/comp);
    
    ratio = ROIAverageBetas_Norm_mean(:,2,iSide)./ROIAverageBetas_Norm_mean(:,1,iSide);
    % ratio = ROIAverageBetas_Norm_mean(:,2,2)./ROIAverageBetas_Norm_mean(:,1,2)
    [ ratio2Plot, level2Plot, fit , error2plot ] =  cal_dbSLvsBetaWeight(ratio,stimulusLevels,baseLevel_dB);
             scatter( level2Plot , ratio2Plot, 50,'filled')
   
    hold on
    fit = polyfit(level2Plot',ratio2Plot,1);
    plot(level2Plot,polyval(fit,level2Plot),'LineWidth',2);
    % correlation = corrcoef([level2Plot_mv' ratio2Plot_mv]);
    % text(25,0.8,sprintf('Correlation = %.2f \n y = %.2fx + %.2f',correlation(2), fit(1), fit(2)))
    % plot(f,f,'k--')
    errorbar(error2plot(1),error2plot(2),error2plot(3))
    % xlabel('Stimulus Sensation Level (dB SL)');
    % ylabel('Average Beta Weight Ratio (ConB / ConA)');
    plot(stimInfo.stimLevel_SL_mv,stimInfo.stimLevel_SL_mv./max(stimInfo.stimLevel_SL_mv),'k--','LineWidth',2)
        title(Info.Sides{iSide},'FontSize',20)
        

    xlim([14 51])
    ylim([0 1.2])
    ax = gca;
ax.FontSize = 15;
    box on
    if iSide == 1
        legend('Data', 'Line of best fit', 'S.T.D. at 50 SL','Normalised SL function','Location','northwest')
        legend('boxoff','FontSize',FontSize)
        ylabel('Average Beta Weight Ratio (sHL / NH)','FontSize',FontSize)
    else
        ax.YTickLabel = []; 
        ax1 = axes('Position',[0.5 0.22 1 1],'Visible','off');
        axes(ax1)
        text(0,0,'Stimulus Sensation Level (dB SL)','HorizontalAlignment','center','FontSize',FontSize)
    end
%     xlabel('Stimulus Sensation Level (dB SL)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%% HERE
conA_color = [0 113 188]./256;
conB_color = [216 82 29]./256;
%% ROI AVERAGE BETA WEIGHT - Both sides
a = cat(4,ROIAverageBetas_Norm(:,:,1,:), ROIAverageBetas_Norm(:,:,2,:));
b = mean(a,4);
e = std(a,0,4) ./sqrt(14);

xlabelnamesMV = stimInfo.stimNames.mv;

figure('Units','centimeters',...
    'Position',[figx bottom width/2 height/2],...
    'color',[1 1 1],...
    'Name','ROI Average Beta Weight')
yyaxis left

errorbar(1:length(b),b(:,1),e(:,1),'color',conA_color)
hold on
errorbar(1:length(b),b(:,2),e(:,2),'color',conB_color)
xlim([1 length(b)])
ylim([0 1])
ax = gca;
ax.YColor = [0 0 0];
ax.XTickLabel = round(xlabelnamesMV(ax.XTick),2);

ylabel('Normalised Beta Weight')
yyaxis right
% plot(stimInfo.stimNames.all,stimInfo.stimLevel_SL/max(stimInfo.stimLevel_SL))
plot(1:length(b),stimInfo.stimLevel_SL_mv/max(stimInfo.stimLevel_SL_mv),...
    'color',[ 0 0 0 ],...
    'linestyle','--')
ylabel('Normalised Sensation Level')
ax = gca;
ax.YColor = [0 0 0];
ylim([0 1])
xlabel('Frequency (kHz)')
% legend('Condition A','Condition B','Sensation Level','Location','northeastoutside')
legend('Condition A','Condition B','Sensation Level','Location','southwest')
legend('boxoff')
% ylim([-150 150])


%%%% TO HERE
%%
ratio = b(:,2)./b(:,1)
% ratio = ROIAverageBetas_Norm_mean(:,2,2) ./ ROIAverageBetas_Norm_mean(:,1,2);
stimulusLevels = stimInfo.stimLevel_SL_mv;
baseLevel_dB = 50;
[ ratio2Plot, level2Plot, fit , error2plot ] =  cal_dbSLvsBetaWeight(ratio,stimulusLevels,baseLevel_dB);

figure('Units','centimeters',...
    'Position',[figx bottom width/2 height/2],...
    'color',[1 1 1],...
    'Name','BOLD vs SL')
scatter( level2Plot , ratio2Plot )
hold on
fit = polyfit(level2Plot',ratio2Plot,1);
plot(level2Plot,polyval(fit,level2Plot));
% correlation = corrcoef([level2Plot_mv' ratio2Plot_mv]);
% text(25,0.8,sprintf('Correlation = %.2f \n y = %.2fx + %.2f',correlation(2), fit(1), fit(2)))
% plot(f,f,'k--')
errorbar(error2plot(1),error2plot(2),error2plot(3))
% xlabel('Stimulus Sensation Level (dB SL)');
% ylabel('Average Beta Weight Ratio (ConB / ConA)');
plot(stimInfo.stimLevel_SL_mv,stimInfo.stimLevel_SL_mv./max(stimInfo.stimLevel_SL_mv),'k--')
legend('Average across ROIs', 'Line of best fit','Location','northwest')
legend('boxoff','FontSize',FontSize)
xlabel('Stimulus Sensation Level (dB SL)'); ylabel('Average Beta Weight Ratio (ConB / ConA)');

%% pRF

for iSub = 1:nSubjects
    
    ROIpCF_conA_left{:,iSub} = subjectData(iSub).data.Left.LeftAC.ConcatenationNH_Smoothed_fwhm2.pRF{1, 2};
    ROIpCF_conB_none_left{:,iSub} = subjectData(iSub).data.Left.LeftAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_None{1, 2};
    ROIpCF_conB_SL_level_left{:,iSub} = subjectData(iSub).data.Left.LeftAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_SL_level{1, 2};
    
    ROIpCF_conA_right{:,iSub} = subjectData(iSub).data.Right.RightAC.ConcatenationNH_Smoothed_fwhm2.pRF{1, 2};
    ROIpCF_conB_none_right{:,iSub}= subjectData(iSub).data.Right.RightAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_None{1, 2};
    ROIpCF_conB_SL_level_right{:,iSub} = subjectData(iSub).data.Right.RightAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_SL_level{1, 2};
    
    % left.ROIpCF_conA_none(:,iSub) = subjectData(iSub).data.Left.LeftAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_None{1, 2};
    % left.ROIpCF_conA_SL_level(:,iSub) = subjectData(iSub).data.Left.LeftAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_SL_level{1, 2};
    
    % right.ROIpCF_conB(:,iSub) = subjectData(iSub).data.Right.RightAC.ConcatenationNH_Smoothed_fwhm2.pRF{1, 2};
    % right.ROIpCF_conB_none(:,iSub) = subjectData(iSub).data.Right.RightAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_None{1, 2};
    % right.ROIpCF_conB_SL_level(:,iSub) = subjectData(iSub).data.Right.RightAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_SL_level{1, 2};
    left.ROI_pRF_A_r2{:,iSub} = subjectData(iSub).data.Left.LeftAC.ConcatenationNH_Smoothed_fwhm2.pRF{1, 1};
    left.ROI_pRF_B_none_r2{:,iSub} = subjectData(iSub).data.Left.LeftAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_None{1, 1};
    left.ROI_pRF_B_SL_r2{:,iSub} = subjectData(iSub).data.Left.LeftAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_SL_level{1, 1};
    
    right.ROI_pRF_A_r2{:,iSub} = subjectData(iSub).data.Right.RightAC.ConcatenationNH_Smoothed_fwhm2.pRF{1, 1};
    right.ROI_pRF_B_none_r2left.ROI_pRF_A_r2{:,iSub}= subjectData(iSub).data.Right.RightAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_None{1, 1};
    right.ROI_pRF_B_SL_r2{:,iSub} = subjectData(iSub).data.Right.RightAC.ConcatenationHLsim_Smoothed_fwhm2.pRF_SL_level{1, 1};
    
end


%% USE
%% histogram
nBins = 8;
for iSide = 1:2
    if iSide == 1
        conAin = ROIpCF_conA_left;
        conBinNone = ROIpCF_conB_none_left;
        conBinSL = ROIpCF_conB_SL_level_left;
        r2 = left.ROI_pRF_A_r2;
        subID = 1:7;
    else
        conAin = ROIpCF_conA_right;
        conBinNone = ROIpCF_conB_none_right;
        conBinSL = ROIpCF_conB_SL_level_right; 
        r2 = right.ROI_pRF_A_r2;
        subID = 8:14;
    end
   binEdges = linspace(0,41.5,9)
   for i = 1:nSubjects
       iSub = subID(i);
       %     histogram(X,edges) sorts X into bins with the bin edges specified by
       % the vector, edges. Each bin includes the left edge, but does not include
       % the right edge, except for the last bin which includes both edges.
       % define percent
       x = r2{:,i};
       Nx = size(x,2);
       kpercent = 25;
       
       % STEP 1 - rank the data
       y = sort(x);
       
       % STEP 2 - find k% (k /100) of the sample size, n.
       k = kpercent/100;
       result = k*Nx;
       
       % STEP 3 - if this is an integer, add 0.5. If it isn't an integer round up.
       [N,D] = rat(k*Nx);
       if isequal(D,1)              % k*Nx is an integer, add 0.5
           result = result+0.5;
       else                           % round up
           result = round(result);
       end
       
       % STEP 4 - Find the number in this position. If your depth ends
       % in 0.5, then take the midpoint between the two numbers.
       [T,R] = strtok(num2str(result),'0.5');
       if strcmp(R,'.5'),
           Qk = mean(y(result-0.5:result+0.5));
       else
           Qk = y(result);
       end

       conAinTemp = conAin{:,i}(r2{:,i}>Qk);
       conBinNoneTemp = conBinNone{:,i}(r2{:,i}>Qk);
       conBinSLTemp = conBinSL{:,i}(r2{:,i}>Qk);
       
       temp = histogram(conAinTemp,binEdges);
       pCFdistrobution_conA(:,iSub) = temp.BinCounts;
       pCFdistrobution_conA_percentage(:,iSub) = (pCFdistrobution_conA(:,i)./(sum(pCFdistrobution_conA(:,i)))).*100;
       BinLimits = temp.BinLimits;
       
       temp = histogram(conBinNoneTemp,binEdges);
       pCFdistrobution_conA_none(:,iSub) = temp.BinCounts;
       pCFdistrobution_conA_none_percentage(:,iSub) = (pCFdistrobution_conA_none(:,i)./(sum(pCFdistrobution_conA_none(:,i)))).*100;
       BinLimits = temp.BinLimits;
       
       temp = histogram(conBinSLTemp,binEdges);
       pCFdistrobution_conA_SL(:,iSub) = temp.BinCounts;
       pCFdistrobution_conA_SL_percentage(:,iSub) = (pCFdistrobution_conA_SL(:,i)./(sum(pCFdistrobution_conA_SL(:,i)))).*100;
       BinLimits = temp.BinLimits;
       BinEdges = temp.BinEdges
       
       
    A = conAinTemp;
    B_none = conBinNoneTemp;
    B_SL = conBinSLTemp;
    
    vectorSumNorm_none(i) = sum((A./mean(A)).*(B_none./mean(B_none)));
    vectorSum_none(i) = sum(A.*B_none);
    vectorSumNorm_SL(i) = sum((A./mean(A)).*(B_SL./mean(B_SL)));
    vectorSum_SL(i) = sum(A.*B_SL);
    
    vectorSumNorm_AvA(i) = sum((A./mean(A)).*(A./mean(A)));
    vectorSum_AvA(i) = sum(A.*A);
    
    
   c_AvA(:,:,i) = corrcoef(A,A);
   c_none(:,:,i)=  corrcoef(A,B_none);
   c_SL(:,:,i) =  corrcoef(A,B_SL);
   end
end
FontSize = 18
bar2plot_mean = [mean(pCFdistrobution_conA_percentage,2),mean(pCFdistrobution_conA_none_percentage,2),mean(pCFdistrobution_conA_SL_percentage,2)];
bar2plot_std = [std(pCFdistrobution_conA_percentage,0,2),std(pCFdistrobution_conA_none_percentage,0,2),std(pCFdistrobution_conA_SL_percentage,0,2)]  ./sqrt(nSubjects*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','centimeters',...
    'Position',[figx bottom width height/1.5],...
    'color',[1 1 1],...
    'Name','Histogram')
h = barwitherr(bar2plot_std,bar2plot_mean);
% set(gca,'XTickLabel',{'Group A','Group B','Group C'})
binEdges_round = round(invNErb(binEdges),1);
set(gca,'XTickLabel',{['> ' mat2str(binEdges_round(2))],[mat2str(binEdges_round(2)) ' > ' mat2str(binEdges_round(3))],...
    [mat2str(binEdges_round(3)) '>' mat2str(binEdges_round(4))],...
    [mat2str(binEdges_round(4)) '>' mat2str(binEdges_round(5))],...
    [mat2str(binEdges_round(5)) '>' mat2str(binEdges_round(6))],...
    [mat2str(binEdges_round(6)) '>' mat2str(binEdges_round(7))],...
    [mat2str(binEdges_round(7)) '>' mat2str(binEdges_round(8))],...
    [mat2str(binEdges_round(8)) '<']})
set(h(1),'FaceColor',conA_color);
set(h(2),'FaceColor',conB_color);
% set(h(3),'FaceColor',[125 46 141]./256);
set(h(3),'FaceColor',[ 46 125 70]./256);
set(h,'BarWidth',1)
ax = gca;
ax.FontSize = 14;
legend('Normal Hearing', 'Hearing Loss Simulation (Standard pRF)',  'Hearing Loss Simulation (Modified pRF)','Location','northeast')
legend('boxoff','FontSize',FontSize)
ylabel('Percentage of voxels','FontSize',FontSize)
xlabel('Binned estimated voxel population centre frequency (kHz)','FontSize',FontSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vectorSum_none_mean = mean(vectorSumNorm_none);
vectorSum_SL_mean = mean(vectorSumNorm_SL);
vectorSum_AvA_mean = mean(vectorSumNorm_AvA);

vectorSum_none_ste = std(vectorSumNorm_none)./sqrt(nSubjects*2);
vectorSum_SL_ste = std(vectorSumNorm_SL)./sqrt(nSubjects*2);
vectorSum_AvA_ste = std(vectorSumNorm_AvA)./sqrt(nSubjects*2);

vectorSum_none_mean = mean(vectorSum_none);
vectorSum_SL_mean = mean(vectorSum_SL);
vectorSum_AvA_mean = mean(vectorSum_AvA);

vectorSum_none_ste = std(vectorSum_none)./sqrt(nSubjects*2);
vectorSum_SL_ste = std(vectorSum_SL)./sqrt(nSubjects*2);
vectorSum_AvA_ste = std(vectorSum_AvA)./sqrt(nSubjects*2);

figure
bar2plot_mean = [vectorSum_AvA_mean vectorSum_none_mean vectorSum_SL_mean];
bar2plot_std = [vectorSum_AvA_ste vectorSum_none_ste vectorSum_SL_ste];
bar(bar2plot_mean)
hold on
errorbar(bar2plot_mean,bar2plot_std,'.')

c_mean = [mean(c_none,3), mean(c_SL,3)]
c_std = [std(c_none,0,3)./sqrt(nSubjects*2) std(c_SL,0,3)./sqrt(nSubjects*2)]

bar2plot_c_mean = [c_mean(2), c_mean(6)]
bar2plot_c_std = [c_std(2), c_std(6)]
figure
bar(bar2plot_c_mean)
hold on
errorbar(bar2plot_c_mean,bar2plot_c_std,'.')

%% NOT USED

% histogram each subject to get nBins x frequency per bin
% divide by number of voxels to get percentage
% mean and ste of each bin
% plot using histogram
% nBins = 8;
% for iSide = 1:2
%     if iSide == 1
%         conAin = ROIpCF_conA_left;
%         conBinNone = ROIpCF_conB_none_left;
%         conBinSL = ROIpCF_conB_SL_level_left;
%     else
%         conAin = ROIpCF_conA_right;
%         conBinNone = ROIpCF_conB_none_right;
%         conBinSL = ROIpCF_conB_SL_level_right;        
%     end
% for iSub = 1:nSubjects
%     temp = histogram(conAin{:,iSub},nBins);
%     pCFdistrobution_conA(:,iSub) = temp.BinCounts;
%     pCFdistrobution_conA_percentage(:,iSub) = (pCFdistrobution_conA(:,iSub)./(sum(pCFdistrobution_conA(:,iSub)))).*100;
%     BinLimits = temp.BinLimits;
%     
%     temp = histogram(conBinNone{:,iSub},nBins);
%     pCFdistrobution_conA_none(:,iSub) = temp.BinCounts;
%     pCFdistrobution_conA_none_percentage(:,iSub) = (pCFdistrobution_conA_none(:,iSub)./(sum(pCFdistrobution_conA_none(:,iSub)))).*100;
%     BinLimits = temp.BinLimits;
%     
%     temp = histogram(conBinSL{:,iSub},nBins);
%     pCFdistrobution_conA_SL(:,iSub) = temp.BinCounts;
%     pCFdistrobution_conA_SL_percentage(:,iSub) = (pCFdistrobution_conA_SL(:,iSub)./(sum(pCFdistrobution_conA_SL(:,iSub)))).*100;
%     BinLimits = temp.BinLimits;
%     
% end
% bar2plot_mean = [mean(pCFdistrobution_conA_percentage,2),mean(pCFdistrobution_conA_none_percentage,2),mean(pCFdistrobution_conA_SL_percentage,2)];
% bar2plot_std = [std(pCFdistrobution_conA_percentage,0,2),std(pCFdistrobution_conA_none_percentage,0,2),std(pCFdistrobution_conA_SL_percentage,0,2)]  ./sqrt(nSubjects);
% figure
% h = barwitherr(bar2plot_std,bar2plot_mean)
% set(h(1),'FaceColor','k');
% % bar(bar2plot_mean)
% % hold on
% % errorbar(bar2plot_mean,bar2plot_std,'.')
% %   h = barwitherr(errY, y);% Plot with errorbars
% %
% %   set(gca,'XTickLabel',{'Group A','Group B','Group C'})
% %   legend('Parameter 1','Parameter 2','Parameter 3','Parameter 4')
% %   ylabel('Y Value')
% %   set(h(1),'FaceColor','k');
% 
% figure
% errorbar(bar2plot_mean,bar2plot_std)
% end
%% now need tocal correlation for each subject then take mean


% histogram(ROIpCF_conA_left{:,iSub},nBins)
% hold on
% histogram(ROIpCF_conA_none_left{:,iSub},nBins)
% histogram(ROIpCF_conA_SL_level_left{:,iSub},nBins)

%% make for loop to concat data and add subject info per voxel
% mean(ROIpCF_conA_left{:,:})
% 
% %% make this for all hemispheres
% for iSide = 1:2
%     if iSide == 1
%         conAin = ROIpCF_conA_left;
%         conBinNone = ROIpCF_conB_none_left;
%         conBinSL = ROIpCF_conB_SL_level_left;
%         subID = 1:7;
%     else
%         conAin = ROIpCF_conA_right;
%         conBinNone = ROIpCF_conB_none_right;
%         conBinSL = ROIpCF_conB_SL_level_right; 
%         subID = 8:14;
%     end
% for iSub = 1:nSubjects
%     i = subID(iSub);
%     
%     A = conAin{:,iSub};
%     B_none = conBinNone{:,iSub};
%     B_SL = conBinSL{:,iSub};
%     
%     vectorSumNorm_none(i) = sum((A./mean(A)).*(B_none./mean(B_none)));
%     vectorSum_none(i) = sum(A.*B_none);
%     vectorSumNorm_SL(i) = sum((A./mean(A)).*(B_SL./mean(B_SL)));
%     vectorSum_SL(i) = sum(A.*B_SL);
%     
%     vectorSumNorm_AvA(i) = sum((A./mean(A)).*(A./mean(A)));
%     vectorSum_AvA(i) = sum(A.*A);
%     
%     
%    [c_AvA(i) p_AvA(i)] = corrcoef(A,A);
%    [c_none(i) p_none(i)] =  corrcoef(A,B_none);
%    [c_SL(i) p_SL(i)] =  corrcoef(A,B_SL);
%    
%  
% end
% end
% 
% vectorSum_none_mean = mean(vectorSumNorm_none);
% vectorSum_SL_mean = mean(vectorSumNorm_SL);
% vectorSum_AvA_mean = mean(vectorSumNorm_AvA);
% 
% vectorSum_none_ste = std(vectorSumNorm_none)./sqrt(nSubjects*2);
% vectorSum_SL_ste = std(vectorSumNorm_SL)./sqrt(nSubjects*2);
% vectorSum_AvA_ste = std(vectorSumNorm_AvA)./sqrt(nSubjects*2);
% 
% vectorSum_none_mean = mean(vectorSum_none);
% vectorSum_SL_mean = mean(vectorSum_SL);
% vectorSum_AvA_mean = mean(vectorSum_AvA);
% 
% vectorSum_none_ste = std(vectorSum_none)./sqrt(nSubjects*2);
% vectorSum_SL_ste = std(vectorSum_SL)./sqrt(nSubjects*2);
% vectorSum_AvA_ste = std(vectorSum_AvA)./sqrt(nSubjects*2);
% 
% figure
% bar2plot_mean = [vectorSum_AvA_mean vectorSum_none_mean vectorSum_SL_mean];
% bar2plot_std = [vectorSum_AvA_ste vectorSum_none_ste vectorSum_SL_ste];
% bar(bar2plot_mean)
% hold on
% errorbar(bar2plot_mean,bar2plot_std,'.')


A_data = [ROIpCF_conA_left, ROIpCF_conA_right];
B_Data_None = [ROIpCF_conB_none_left, ROIpCF_conB_none_right];
B_Data_SL = [ROIpCF_conB_SL_level_left, ROIpCF_conB_SL_level_right];

figure
for i = 1:nSubjects*2
    subplot(1,2,1)
    scatter(A_data{i},B_Data_None{i})
    hold on
    
    fit_none(i,:) = polyfit(A_data{i},B_Data_None{i},1);
    plot(1:45,polyval(fit_none(i,:),1:45));
    subplot(1,2,2)
    scatter(A_data{i},B_Data_SL{i})
    hold on
    fit_SL(i,:) = polyfit(A_data{i},B_Data_SL{i},1);
    plot(1:45,polyval(fit_SL(i,:),1:45));
end



% % correlation
% index = ~isnan(conA_data{2});
% % vectorSum = sum(conA_data{2}(index).*conB_data{2}(index));
% A = conA_data{2}(index);
% B = conB_data{2}(index);
% data.vectorSumNorm = sum((A./mean(A)).*(B./mean(B)));
% data.vectorSum = sum(A.*B);

%
%
% plot(mean(ROIAverageTuningCurves_Norm(:,:,2,:),4))
% max([max(max(ROIAverageTuningCurves_ConA(:,:,1))), max(max(ROIAverageTuningCurves_ConB(:,:,1)))])
%
% figure
% % ROIAverageTuningCurves_Norm

% g = gramm('x',  ,'y',data1.b,'color',data1.Q);
% % g.geom_point('alpha',0.05)
% g.facet_grid([],data1.Q);
% % g.facet_wrap(data.Q);
% g.geom_point()
%
% % g.stat_smooth();
% g.stat_cornerhist();
% % g.stat_glm();
% % g.stat_bin();
% g.axe_property('DataAspectRatio',[1 1 1]);
% % Set appropriate names for legends
% g.set_names('column','Quantile','x','Condition A','y','Condition B','color','r2');
% %%%
% % Set figure title
% g.set_title('Estimated pCF');
% figure
% g.draw()

% pRF
% voxel pCF comparisions (scatter plot)
% voxel pCF distribution

%% pre process data
% normalise
%   divide by maximium value in condition A (normal hearing)

%% Plot data

end