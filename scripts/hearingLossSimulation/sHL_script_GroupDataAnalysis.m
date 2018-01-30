function groupData = sHL_script_GroupDataAnalysis
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


% ROIAverageTuningCurves_Norm(:,:,:,5) = nan;

xlabel = stimInfo.stimNames.bin;
for iSide = 1:2
    figure
    
minBeta = min(min(min(min(ROIAverageTuningCurves_Norm(:,:,:,iSide,:)))));
maxBeta = max(max(max(max(ROIAverageTuningCurves_Norm(:,:,:,iSide,:)))));
for i = 1:8
subplot(2,4,i)
meanBetasA = nanmean(ROIAverageTuningCurves_Norm(i,:,1,iSide,:),5);
meanBetasB = nanmean(ROIAverageTuningCurves_Norm(i,:,2,iSide,:),5);

stdBetasA = nanstd(ROIAverageTuningCurves_Norm(i,:,1,iSide,:),0,5) ./sqrt(nSubjects);
stdBetasB = nanstd(ROIAverageTuningCurves_Norm(i,:,2,iSide,:),0,5) ./sqrt(nSubjects);
% 
% minBeta = min([meanBetasA, meanBetasB]) + min([stdBetasA, stdBetasB]);
% maxBeta = max([meanBetasA, meanBetasB]) + max([stdBetasA, stdBetasB]);
% stderror = std( data ) / sqrt( length( data ))
% plot(nanmean(ROIAverageTuningCurves_Norm(i,:,1,:),4))
errorbar(meanBetasA,stdBetasA)
hold on
% plot(nanmean(ROIAverageTuningCurves_Norm(i,:,2,:),4))
errorbar(meanBetasB,stdBetasB)

xlim([1 8])
ylim([minBeta maxBeta])
    ax = gca;
ax.XTickLabel = round(xlabel(ax.XTick),2);
end
end

% GLM ROI average


plot(mean(ROIAverageTuningCurves_Norm(:,:,2,:),4))
max([max(max(ROIAverageTuningCurves_ConA(:,:,1))), max(max(ROIAverageTuningCurves_ConB(:,:,1)))])

figure
% ROIAverageTuningCurves_Norm

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