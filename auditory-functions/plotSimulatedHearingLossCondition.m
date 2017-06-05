function plotSimulatedHearingLossCondition(frq)
frq = logspace(log10(0.02),log10(20),500);
NLevel = 25;
CriticalRatioFFT = getCriticalRatioPerERB(frq*1000);
thresholdBaselineFFT = NLevel * ones(size(frq));
threshold = 75;
thresholdHearingLossFFT = funSimulateHearingLoss(frq);
thresholdHearingLossFFT = min(thresholdHearingLossFFT,threshold);
levelFFT = max(thresholdHearingLossFFT,thresholdBaselineFFT) - CriticalRatioFFT;

allFrequencies = funInvNErb(linspace(funNErb(0.1),funNErb(8),32));

allFrequenciesLevel = repmat(75,length(allFrequencies),1);

%     plot(allFrequencies,allFrequenciesLevel)
allUpperF = funInvNErb(funNErb(allFrequencies) + 0.5);
allLowerF = funInvNErb(funNErb(allFrequencies) - 0.5);
%
stimuli = zeros(size(frq));
%     stimuli(frq(allFrequencies)) = 75;
stimIm = zeros(length(allFrequencies),length(frq));
for i = 1:length(allFrequencies)
    %         stimuli(frq<allUpperF(i) & frq>allLowerF(i)) = 75;
    [val index]  =  min(abs(frq-allFrequencies(i)));
    index2  = frq<allUpperF(i) & frq>allLowerF(i);
    stimuli(index) = 75;
    stimIm (i,index2) = 75;
 stimImgramm{i} = stimIm(i,:);
 stimIDgramm{i} = repmat(i,size(stimIm(i,:)))
end

% set up variables from gramm
x = repmat(frq,1,4);
y = [thresholdHearingLossFFT CriticalRatioFFT levelFFT stimuli];
names = [repmat({'tsHL (dBSPL)'},length(frq),1); repmat({'CR (dBSPL)'},length(frq),1); repmat({'sHL (dBSPL)'},length(frq),1) ; repmat({'Stimuli dBSPL'},length(frq),1)];
type = [repmat({'Background Noise'},length(frq)*3,1); repmat({'stimuli'},length(frq),1)];
figure
g(1,1) = gramm('x',x,'y',y,'color',names,'marker',type);
g(1,1).geom_line()
% Set appropriate names for legends
g(1,1).set_names('x','Frequency (kHz)','y','Level (dB)','color','Level:')
%Set figure title
g(1,1).set_title('Hearing loss simulation')

g.draw()


% [xALimMin xALimMinIn] = min(abs(frq-0.05));
% [xALimMax xALimMaxIn] = min(abs(frq-10));

%%%%% PLOT STIMULI %%%%%
figure('Color',[1 1 1]);
for ii = 1:length(allFrequencies)
    plotColor = [allFrequencies(ii)/max(allFrequencies) 0.4 0.6];
%     semilogx(frq,stimIm(ii,:), 'color',plotColor)
    plot(frq,stimIm(ii,:), 'color',plotColor)    
% xlim([frq(xALimMinIn) frq(xALimMaxIn)])
    xlim([0.05 12])
    set(gca,'XScale','log','XTick',[0.1 1 10]);
    hold on
    fill(frq,stimIm(ii,:),plotColor)
%     l = legend([l num2str(ii)])
end


figure('Color',[1 1 1]);
g = gramm('x',frq,'y',stimImgramm,'color',stimIDgramm);
g.geom_line()
% Set appropriate names for legends
g.set_names('x','Frequency (kHz)','y','Level (dB)','color','Stimulus ID:')
g.axe_property('XScale','log');
g.axe_property('YLim',[0 80]);
g.axe_property('XLim',[0.05 10]);
g.draw()

% xlim([frq(xALimMinIn) frq(xALimMaxIn)])

%%%%% PLOT Con A %%%%%
figure('Color',[1 1 1]);
semilogx(frq,thresholdBaselineFFT,'c', 'LineWidth', 2);
hold on
semilogx(frq,CriticalRatioFFT,'r', 'LineWidth', 2)
semilogx(frq,thresholdBaselineFFT - CriticalRatioFFT, 'b', 'LineWidth', 3)
for ii = 1:length(allFrequencies)
    semilogx(allFrequencies(ii),75,'ro','MarkerFaceColor','r')
    hold on
end
xlim([min(frq) max(frq)])

%%%%% PLOT Con B %%%%%
figure('Color',[1 1 1]);
semilogx(frq,thresholdBaselineFFT,'c', 'LineWidth', 2)
hold on
semilogx(frq,thresholdHearingLossFFT,'g', 'LineWidth', 2)
semilogx(frq,CriticalRatioFFT,'r', 'LineWidth', 2)
semilogx(frq,max(thresholdHearingLossFFT,thresholdBaselineFFT) - CriticalRatioFFT, 'b', 'LineWidth', 3)
for ii = 1:length(allFrequencies)
    semilogx(allFrequencies(ii),75,'ro','MarkerFaceColor','r')
    hold on
end
xlim([min(frq) max(frq)])



%%%%%
figure('Color',[1 1 1]);
%         plot(frq,ee)
hold on
plot(frq,thresholdBaselineFFT,'g')
plot(frq,thresholdHearingLossFFT,'b')
plot(frq,CriticalRatioFFT,'r')
plot(frq,max(thresholdHearingLossFFT,thresholdBaselineFFT) - CriticalRatioFFT, 'm', 'LineWidth', 2)

for ii = 1:length(allFrequencies)
    plot(allFrequencies(ii),75,'ro')
    hold on
end

xlim([min(frq) max(frq)])
xlabel('Frequency (kHz)')
ylabel('dB SPL')
legend( 'Baseline threshold',...
    'Hearing Loss threshold',...
    'Critical Ratio',...
    'max(thresholdHearingLossFFT,thresholdBaselineFFT) - CriticalRatioFFT',...
    'Location','best')
title('Gain applied to EE noise (Normalised to 1Vrms in 1 ERB) as a function of Frequency')
set(gcf,'color','w');