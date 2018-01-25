function data = script_pRFROIAnalysis(conA_data,conB_data,dataNames,analysisName)

%% could put all data for all analysis in here and compare using subplot/facets

% correlation vs cona frequency??

%% put data into gramm struct
% remove nans

index = ~isnan(conA_data{2});
% vectorSum = sum(conA_data{2}(index).*conB_data{2}(index));
A = conA_data{2}(index);
B = conB_data{2}(index);
data.vectorSumNorm = sum((A./mean(A)).*(B./mean(B)));
data.vectorSum = sum(A.*B);
%% repmat dataNames to match each data cell
data1.names = repmat(dataNames,size(conA_data));
% figure
% scatter(conA_data{2},conB_data{2})
% xlim = ([0 45]);
% ylim = ([0 45]);
% figure
% histogram(conA_data{1})

% histogram
data1.a = conA_data{2};
data1.b = conB_data{2};
% data.r2 = conA_data{1} + conB_data{1};

data1.r2 = conA_data{1};

h = histogram(conA_data{1},4);

Y = prctile(conA_data{1},[25 50 75]);
Q = quantile(conA_data{1},[0.25 0.50 0.75]);
Q = [0 Q 1];
data1.Q = nan(size(data1.r2));

for i = 1:length(Q)-1
data1.Q(data1.r2>=Q(i) & data1.r2<=Q(i+1)) = i;
end


% g = gramm('x',data.a,'y',data.b,'lightness',data.Q);

g = gramm('x',data1.a,'y',data1.b,'color',data1.Q);
% g.geom_point('alpha',0.05)
g.facet_grid([],data1.Q);
% g.facet_wrap(data.Q);
g.geom_point()

% g.stat_smooth();
g.stat_cornerhist();
% g.stat_glm();
% g.stat_bin();
g.axe_property('DataAspectRatio',[1 1 1]);
% Set appropriate names for legends
g.set_names('column','Quantile','x','Condition A','y','Condition B','color','r2');
%%%
% Set figure title
g.set_title('Estimated pCF');
figure
g.draw()



data2.pCF = [conA_data{2},conB_data{2}];
data2.r2 = [conA_data{1},conB_data{1}];
data2.r2_a = [conA_data{1},conA_data{1}];
data2.r2_diff = [conA_data{1} - conB_data{1},conA_data{1} - conB_data{1}];
data2.condition = [repmat({'Condition A'}, size(conA_data{2})), repmat({'Condition B'}, size(conB_data{2}))];
Q = quantile(data2.r2_a,[0.25 0.50 0.75]);
Q = [0 Q 1];
data2.Q = nan(size(data2.r2));

for i = 1:length(Q)-1
data2.Q(data2.r2>=Q(i) & data2.r2<=Q(i+1)) = i;
end

clear g
g = gramm('x',data2.pCF,'color',data2.condition);
% g.geom_point('alpha',0.05)
% g.facet_grid([],data.Q);
% g.geom_point()
g.stat_bin('nbins',8,'geom','overlaid_bar');
% g.axe_property('DataAspectRatio',[1 1 1]);
% Set appropriate names for legends
g.set_names('column','Quantile','x','Frequency (nERBs)','y','Number of Voxels','color','r2');
%%%
% Set figure title
g.set_title('Estimated pCF');
figure
g.draw()
end
% create function or add to: voxel comparisions;
% pCF scatter, pCF distribution, pCF correlation
% r = corr2(A,B)


% 
% max(data(iGroup).r2{i}(r2Index))
% for i = 1:2
% figure
% 
% VoxelStruct(i).conA = data(1).pCFest{i}(data(2).pCFest{1}>3 & data(2).pCFest{1}<33);
% VoxelStruct(i).conB = data(2).pCFest{1}(data(2).pCFest{1}>3 & data(2).pCFest{1}<33);
% 
% g(i) = gramm('x',VoxelStruct(i).conA,'y',VoxelStruct(i).conB);
% % g.geom_point('alpha',0.05)
% g(i).geom_point()
% g(i).stat_glm();
% 
% % g(i).stat_bin();
% g(i).draw()
% end
% 
% % nhIndex = data(2).pCFest{1}>3 & data(2).pCFest{1}<33;
% nhIndex = length(data(2).pCFest{1});
% for i = 1:2
% figure
% VoxelStruct(i).conB = data(1).pCFest{i}(nhIndex & r2Index);
% VoxelStruct(i).conA = data(2).pCFest{1}(nhIndex & r2Index);
% g(i) = gramm('x',VoxelStruct(i).conA,'y',VoxelStruct(i).conB);
% % g.geom_point('alpha',0.05)
% g(i).geom_point()
% g(i).stat_glm();
%  g(i).set_names('x','Condition A - Frequency (ERB)','y','Condition B - Frequency (ERB)')
% % g(i).stat_bin();
% g(i).draw()
% end
% 
% unpacked.pCF = [data(1).pCFest{1}(r2Index) data(1).pCFest{2}(r2Index) data(2).pCFest{1}(r2Index)];
% % unpacked.con = [repmat(1,1,length(data(1).pCFest{1})) repmat(2,1,length(data(1).pCFest{1})) repmat(3,1,length(data(1).pCFest{1}))];
% unpacked.con = [repmat({'Condition B - pRF'},1,sum(r2Index)) repmat({'Condition B - pRF_mod'},1,sum(r2Index)) repmat({'Condition A'},1,sum(r2Index))];
% 
% 
% figure
% hist  = gramm('x',unpacked.pCF,'color',unpacked.con);
% hist.stat_bin('nbins',8);
% hist.set_names('x','Frequency (ERB)','y','Number of Voxels','color','Condition:')
% % hist.geom_label()
% hist.draw()
% 
% hist = [];
% figure
% hist  = gramm('x',unpacked.pCF,'color',unpacked.con);
% hist.set_names('x','Frequency (ERB)','y','Frequency','color','Condition:')
% hist.stat_density();
% hist.draw()
% 
% % violin pCF vs pTW?
% 
% %%%%
% % make stimulus plot into a gramm object
% figure
% x = 0:40;
% mu = pCFest{1}(2000);
% sigma = pTWest{1}(2000);
% rfModel = exp(-(((x-mu).^2)/(2*(sigma^2))));