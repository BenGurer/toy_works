function plot_compareConditions_ROIav(conA,conB,conID,figureName)
x = 1:(length(conA));
yLimMax = max(max([conA,conB]));
yLimMin = min(min([conA,conB,conA - conB]));

figure('Name',[figureName, '- Average Beta Weights'],'color',[1 1 1])

plot(x,conA,'color',[0 113 188]/255)
hold on
plot(x,conB,'color',[216 82 24]/255)

ylim([yLimMin yLimMax])
xlim([min(x) max(x)])

xlabel('Condition (kHz)')
ylabel('Beta weight')
title('ROI average Beta weights')
plot(x,conA - conB,'--r')
plot(x,zeros(1,length(conA)),'--k')
legend('Condition A','Condition B','A - B')
ax = gca;
ax.XTickLabel = round(conID(ax.XTick),2);
end