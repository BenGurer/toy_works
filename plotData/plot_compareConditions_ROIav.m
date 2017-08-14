function data = plot_compareConditions_ROIav(conA,conB,nCons)

figure('color',[1 1 1])
yLimMax = max(max([conA,conB]));
yLimMin = min(min([conA,conB]));

plot(1:length(conA),conA,'color',[0 113 188]/255)
hold on
plot(1:length(conA),conB,'color',[216 82 24]/255)
ylim([yLimMin yLimMax])
xlim([1 length(conA)])
xlabel('Beta')
ylabel('Beta weight')
title('ROI average Beta weights')
plot(1:length(conA),conA - conB,'--r')
plot(1:length(conA),zeros(1,length(conA)),'--k')
legend('Condition A','Condition B','A - B')
end