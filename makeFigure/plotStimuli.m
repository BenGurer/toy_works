function plotStimuli
stimTR = 7500;
TR = 7500;
[params,stimulus] = corticalMagnification([],[],stimTR,TR);
scanner.frequency = 1;
scanner.bandwidth = 100;
scanner.duration = 1250;

% stimulus(1).frequency
% stimulus(1).duration
% stimulus(1).number

allFrequencies = funInvNErb(linspace(funNErb(params.lowFrequency),funNErb(params.highFrequency),params.nFrequencies));
lowCuttingFrequencies = funInvNErb(funNErb(allFrequencies)-params.bandwidthERB/2);
highCuttingFrequencies = funInvNErb(funNErb(allFrequencies)+params.bandwidthERB/2);
allFrequencies = (lowCuttingFrequencies+highCuttingFrequencies)/2;
allBandwidths = (highCuttingFrequencies-lowCuttingFrequencies);
lowCuttingFrequencies = round(lowCuttingFrequencies*1000);
highCuttingFrequencies = round(highCuttingFrequencies*1000);

f = 50:10:10000;
silence = -10;
% stimuli = nan(length(f),length(allFrequencies));
stimuli = repmat(silence,length(f),length(allFrequencies));
for i = 1:length(allFrequencies)
%     indexLower = f(f>=lowCuttingFrequencies(i) & f<=highCuttingFrequencies(i));    
    index= f>=lowCuttingFrequencies(i) & f<=highCuttingFrequencies(i);
stimuli(index,i) = i;

stimuliFrequencies{i} = f(index);


end
% figure;surf(f,1:length(allFrequencies),stimuli)
% 
% figure;surf(stimuli)
% figure;
% for i = 1:size(stimuli,2)
% imagesc(f,i,stimuli(:,i))
% hold on
% end
% figure;semilogx(f,stimuli)
% plotSignal(stimuli)
% plotStimulusImage = x = time;y = frequency;
%
% x = time;
% y = frequency;
% stimImageY = zeros(length(lowCuttingFrequencies) + length(highCuttingFrequencies),1);
% stimImageY = zeros(length(allFrequencies),1);
c = 0
for i = 1:6
    id = stimulus(i).number;       
    tempY = stimImageY;
%     tempY(id:id+1) = 1;
%     stimImage{i} = [ones(length(allFrequencies),scanner.duration) zeros(length(allFrequencies),stimulus(i).duration(1)-scanner.duration)];
%     for ii = 1:length(stimulus(i).duration)-1
%         if isnan(stimulus(i).frequency(ii+1))
%             stimImage{i} = [stimImage{i}, repmat(zeros(length(allFrequencies),1),1,stimulus(i).duration(ii+1))];
%         else
%             stimImage{i} = [stimImage{i}, repmat(tempY,1,stimulus(i).duration(ii+1))];
%         end        
%     end
    
        stimImage{i} = [ones(length(f),scanner.duration) repmat(silence,length(f),stimulus(i).duration(1)-scanner.duration)];
        
    for ii = 1:length(stimulus(i).duration)-1
        c = c + scanner.duration;
        if isnan(stimulus(i).frequency(ii+1))
            stimImage{i} = [stimImage{i}, repmat(silence,length(f),stimulus(i).duration(ii+1))];
            
    stimIm{i}{ii} = 0;
    stimCon{i}{ii} = id;
    stimTime{i}{ii} = 0;
    stimDur{i}{ii} = c + 50;
        else
            stimImage{i} = [stimImage{i}, repmat(stimuli(:,id),1,stimulus(i).duration(ii+1))];
            
    stimIm{i}{ii} = stimuliFrequencies{id};
    stimCon{i}{ii} = repmat(id,1,length(stimuliFrequencies{id}));
    stimTime{i}{ii} = repmat(stimulus(i).duration(ii)+c,1,length(stimuliFrequencies{id}));
    stimDur{i}{ii} = repmat(stimulus(i).duration(ii),1,length(stimuliFrequencies{id}));
        end        
    end
    c = c +7500;
    stimImage{i} = stimImage{i}*100;
    

   
end

figure
g = gramm('x',stimTime{1},'y',stimIm{1});
g.geom_point()
g.draw()

figure;imagesc(stimImage{1})

figure;imagesc(1:length(stimImage{1})*6, f,[stimImage{1} stimImage{2} stimImage{3} stimImage{4} stimImage{5} stimImage{6}])
YScale = 'log';
colormap jet
end
  % ***** plotSignal *****
  function [hCursorT, hCursorF] = plotSignal(signal)
    
    signal = signal(1,:);
    time = sampleDuration*(0:length(signal)-1)/1000;
    hold(hTimeseries,'off');
    plot(hTimeseries,time,signal,'k')
    hold(hTimeseries,'on');
    plot(hTimeseries,[0 TDTcycle/1000],[maxVoltage maxVoltage],'r--');
    plot(hTimeseries,[0 TDTcycle/1000],-1*[maxVoltage maxVoltage],'r--');
    hCursorT = plot(hTimeseries,[0 0],get(hTimeseries,'Ylim'),'r');
    if nStimTRs>1
      plot(hTimeseries,repmat((1:nStimTRs-1)*stimTR/1000,2,1), repmat(get(hTimeseries,'Ylim')',1,nStimTRs-1),'r:');
    end
    set(hTimeseries,'XLim',[0 TDTcycle/1000])
    set(get(hTimeseries,'XLabel'),'String','Time(s)','FontName','Arial','FontSize',10)
    set(get(hTimeseries,'YLabel'),'String','Amplitude','FontName','Arial','FontSize',10)
    set(hTimeseries,'XLim',[0 TDTcycle/1000]);
    title(hTimeseries,'Cursor position is approximate !');

    [specg,frq,t] = spectrogramFunction(signal,round(100/sampleDuration),round(80/sampleDuration),round(50/sampleDuration),1000/sampleDuration);
    hold(hSpectrogram,'off');
    surf(hSpectrogram,t,frq,10*log10(abs(specg)),'EdgeColor','none');
    hold(hSpectrogram,'on');
    grid(hSpectrogram,'off');
    box(hSpectrogram,'on');
    view(hSpectrogram,0,90); 
    axis(hSpectrogram,'tight'); 
    colormap(hSpectrogram,jet); 
    hCursorF = plot(hSpectrogram,[0 0],get(hSpectrogram,'Ylim'),'k');
    set(hSpectrogram,'XLim',[0 TDTcycle/1000]);
    if nStimTRs>1
      plot(hSpectrogram,repmat((1:nStimTRs-1)*stimTR/1000,2,1), repmat(get(hSpectrogram,'Ylim')',1,nStimTRs-1),'k:');
    end
    set(get(hSpectrogram,'XLabel'),'String','Time (s)','FontName','Arial','FontSize',10)
    set(get(hSpectrogram,'YLabel'),'String','Frequency (Hz)','FontName','Arial','FontSize',10)

  end