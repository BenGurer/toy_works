function data = plot_compareConditions_pTW(conA,conB,nfperRow,stimFreqs)

figure('color',[1 1 1])
subIndex = [size(conA,1)/(size(conA,1)/nfperRow) size(conA,1)/nfperRow];
yLimMax = max(max([conA,conB]));
yLimMin = min(min([conA,conB]));
for i = 1:length(conA)
    
    subplot(subIndex(1),subIndex(2),i)    
    plot(stimFreqs,conA(i,:),'color',[0 113 188]/255)
    hold on
    plot(stimFreqs,conB(i,:),'color',[216 82 24]/255)
    ylim([yLimMin yLimMax])
    xlim([min(stimFreqs) max(stimFreqs)])
%     text(stimFreqs(i),conA(i,i),sprintf('%.2f kHz',stimFreqs(i)))
    title([num2str(stimFreqs(i)) ' kHz'])
    plot(stimFreqs,zeros(1,length(conA)),'--k')
    if i == 1
    legend('conA','conB')
    end
    
    set(gca, 'XScale', 'log')
end

end