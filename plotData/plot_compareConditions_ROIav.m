function data = plot_compareConditions_ROIav(conA,conB,conID)

figure('color',[1 1 1])
yLimMax = max(max([conA,conB]));
yLimMin = min(min([conA,conB,conA - conB]));

plot(conID,conA,'color',[0 113 188]/255)
hold on
plot(conID,conB,'color',[216 82 24]/255)
ylim([yLimMin yLimMax])
% xlim tight
xlabel('Condition (kHz)')
ylabel('Beta weight')
title('ROI average Beta weights')
plot(conID,conA - conB,'--r')
plot(conID,zeros(1,length(conA)),'--k')
legend('Condition A','Condition B','A - B')
end