function plot_compareConditions_pTW(conA,conB,nRows,stimIDs,figureName)
figure('Name',[figureName, '-pTW'],'color',[1 1 1])
subIndex = round ([size(conA,1)/(size(conA,1)/nRows) size(conA,1)/nRows]);
yLimMax = max(max([conA,conB]));
yLimMin = min(min([conA,conB]));
x = 1:(size(conA,2));
for i = 1:length(conA)
    subplot(subIndex(1),subIndex(2),i)
    plot(x,conA(i,:),'color',[0 113 188]/255)
    hold on
    plot(x,conB(i,:),'color',[216 82 24]/255)
    ylim([yLimMin yLimMax])
    xlim([min(x), max(x)])
    title([num2str(round(stimIDs(i),2)) ' kHz'])
    xlabel('Frequency (kHz)')
    ylabel('Beta weight')
    plot(x,zeros(1,length(conA)),'--k')
    if i == 1
        legend('conA','conB')
    end
    ax = gca;
ax.XTickLabel = round(stimIDs(ax.XTick),2);
end


% figure('color',[1 1 1])
% contour(x,x,conA)
% ax = gca;
% ax.XTickLabel = round(stimIDs,2);
% ax.YTickLabel = round(stimIDs,2);
% ylim([min(x), max(x)])
% xlim([min(x), max(x)])
% [U,V] = gradient(conA,0.1,0.1);
% hold on
% quiver(x,x,U,V)
% hold off
% 
% figure('color',[1 1 1])
% contour(x,x,conB)
% ax = gca;
% ax.XTickLabel = round(stimIDs,2);
% ax.YTickLabel = round(stimIDs,2);
% ylim([min(x), max(x)])
% xlim([min(x), max(x)])
% [U,V] = gradient(conA,0.1,0.1);
% hold on
% quiver(x,x,U,V)
% hold off
end