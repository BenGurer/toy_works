function data = plot_compareConditions_pTW(conA,conB,nfperRow)

figure('color',[1 1 1])
subIndex = [size(conA,1)/(size(conA,1)/nfperRow) size(conA,1)/nfperRow];
yLimMax = max(max([conA,conB]));
yLimMin = min(min([conA,conB]));
for i = 1:length(conA)
    
    subplot(subIndex(1),subIndex(2),i)    
    plot(1:length(conA),conA(i,:),'color',[0 113 188]/255)
    hold on
    plot(1:length(conA),conB(i,:),'color',[216 82 24]/255)
    ylim([yLimMin yLimMax])
    xlim([1 length(conA)])
    
    plot(1:length(conA),zeros(1,length(conA)),'--k')
end

end